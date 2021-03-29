function options = fmrwhy_workflow_qc(settings_fn)
% A custom workflow that does anatomical and functional data quality control for all subjects in a BIDS directory.
% Steps include anatomical-to-functional processing, basic functional time series preprocessing,
% generating a range of QC metrics and images, and compiling an HTML-report per subject.
% 
% :param settings_fn: User preference settings for the specific analysis and dataset, see 
% :type settings_fn: M-file
% :returns: ``options`` - the updated structure with parameter settings pertaining to the analysis and dataset
% 
% :Outputs: A directory with all assets necessary for rendering the HTML report, located in ``derivatives/fmrwhy-qc/sub-XXX/report_[yyyymmddhhmmss]``
%
% :Example:
%
% >>> settings_fn = '/fMRwhy/fmrwhy/settings/fmrwhy_settings_template.m';
% >>> options = fmrwhy_workflow_qc(settings_fn)
% 
% .. seealso::
%   - A thorough description of the :ref:`quality_reporting` pipeline
%   - A sample report can be `viewed here`_.
% 
% .. _viewed here: https://jsheunis.github.io/fmrwhy_sample_QCreport.html

options = fmrwhy_defaults();

% -------
% SETUP STEP A -- Check dependencies, Matlab path, etc
% -------
options = fmrwhy_util_checkDependencies(options);

% -------
% SETUP STEP B -- Load settings, defaults, filenames and parameters
% -------
% Run settings file ==> populates study/data-specific fields in the options structure, including BIDS variables
run(settings_fn);

% Setup fmrwhy derivative directories on workflow level
options = fmrwhy_bids_setupQcDerivDirs(options.bids_dir, options);

% Validate settings
options = fmrwhy_settings_validate(options)

% Load the subjects
subs = options.subjects_output;

% -------
% QC Pipeline -- for each subject
% -------

for s = 1:numel(subs)
    sub = subs{s};
    % Setup fmrwhy derivatives directories on subject level (this copies data from the main bids_dir)
    options = fmrwhy_bids_setupQcSubDirs(options.bids_dir, sub, options);
 
    % Update workflow options with subject anatomical derivative filenames
    options = fmrwhy_bids_getAnatDerivs(options.bids_dir, sub, options);
 
    % -------
    % STEP 0.2 -- Create functional template
    % -------
    % Create, if it does not exist
    [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', sub, 'task', options.template_task, 'run', options.template_run, 'space', 'individual', 'ext', '_bold.nii')
    template_fn = fullfile(options.preproc_dir, filepath, filename);
    if ~exist(template_fn, 'file')
        disp(['Template functional image does not exist yet. Creating now: ' template_fn]);
        [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', sub, 'ses', options.template_session, 'task', options.template_task, 'run', options.template_run, 'echo', options.template_echo, 'ext', '_bold.nii')
        functional_fn = fullfile(options.preproc_dir, filepath, filename);
        fmrwhy_util_saveNiftiFrom4D(functional_fn, template_fn, 1)
    else
        disp(['Template functional image exists: ' template_fn]);
    end
    options.template_fn = template_fn;
 
 
    % -------
    % STEP 1 -- Structural-functional preprocessing: fmrwhy_preproc_structFunc.m
    % -------
    % Loop through all standard structFunc output filenames and see if these files exist
    struct_func_out_fns = [{options.coregest_anatomical_fn} options.probseg_fns options.transform_fns options.rall_fns options.mask_fns];
    run_structFunc = 0;
    for i = 1:numel(struct_func_out_fns)
        if ~exist(struct_func_out_fns{i}, 'file')
            disp(['Structural-funcional preprocessing output file does not exist yet: ' struct_func_out_fns{i}]);
            run_structFunc = 1;
        end
    end
    % If some of the files do not exist, run the fmrwhy_preproc_structFunc processing pipeline
    if run_structFunc
        disp('Running complete structural-funcional preprocessing pipeline')
        fmrwhy_preproc_structFunc(options);
        disp('Complete!')
        disp('---')
    else
        disp('Structural-funcional preprocessing already completed.')
        disp('---')
    end
 
    % -------
    % STEP 2 -- Anatomical localiser: fmrwhy_preproc_anatLocaliser.m
    % TODO: add more checks to see if this was already done, and add logic to decide what to do
    % -------
    if options.map_rois == 1
        anatLocaliser_fns = {};
        for i = 1:numel(options.tasks)
            % Ignore the 'rest' task (assume there is no task ROI for this; have to change in future if RSnetworks available to be normalised or something)
            if strcmp(options.tasks{i}, 'rest') ~= 1
                % Loop through all ROIs for the particular task
                for j = 1:numel(options.roi.(options.tasks{i}).orig_fn)
                    rroi_fn = options.roi.(options.tasks{i}).rroi_fn{j};
                    [filename, filepath] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', options.roi.(options.tasks{i}).desc{j}, 'ext', '_roi_montage.png');
                    montage_fn = fullfile(options.qc_dir, filepath, filename);
                    anatLocaliser_fns = [anatLocaliser_fns {rroi_fn, montage_fn}];
                end
            end
        end
        run_anatLocaliser = 0;
        for i = 1:numel(anatLocaliser_fns)
            if ~exist(anatLocaliser_fns{i}, 'file')
                disp(['Anatomical localiser output file does not exist yet: ' anatLocaliser_fns{i}]);
                run_anatLocaliser = 1;
            end
        end
        % If some of the files do not exist, run the fmrwhy_preproc_anatLocaliser processing pipeline
        if run_anatLocaliser
            fmrwhy_preproc_anatLocaliser(sub, options)
            disp('Complete!')
            disp('---')
        else
            disp('Anatomical localiser processing already completed.')
            disp('---')
        end
    end
 
    % -------
    % PER SESSION, TASK and RUN
    % -------
 
    % Loop through sessions, tasks, runs
    % options.sessions;
    % options.tasks;
    % options.runs;
    % options.has_sessions;
    % options.has_runs;

    % Find sessions for the current subject
    sub_sessions = bids.query(options.bids_dataset,'sessions','sub', sub);

    % If there are sessions, loop through them:
    if ~isempty(sub_sessions)
        for ses = 1:numel(sub_sessions)
            session = sub_sessions{ses};
            tasks = bids.query(options.bids_dataset,'tasks','sub', sub, 'ses', ses);
            for t = 1:numel(tasks)
                task = tasks{t};
                % Get runs for current task
                runs = bids.query(options.bids_dataset,'runs','sub', sub, 'ses', ses, 'task', task);
                % If there are runs, loop through them
                if ~isempty(runs)
                    for r = 1:numel(runs)
                        current_run = runs{r};
                        % STEP 1 -- Basic functional preprocessing: fmrwhy_bids_preprocFunc.m
                        fmrwhy_bids_preprocFunc(options.bids_dir, sub, task, options, 'ses', session, 'run', current_run);
                        % STEP 2 -- Quality control pipeline: fmrwhy_bids_qcRun.m
                        fmrwhy_bids_qcRun(options.bids_dir, sub, task, options, 'ses', session, 'run', current_run, 'echo', options.template_echo);
                    end
                else % If there are NOT runs:
                    % STEP 1 -- Basic functional preprocessing: fmrwhy_bids_preprocFunc.m
                    fmrwhy_bids_preprocFunc(options.bids_dir, sub, task, options, 'ses', session);
                    % STEP 2 -- Quality control pipeline: fmrwhy_bids_qcRun.m
                    fmrwhy_bids_qcRun(options.bids_dir, sub, task, options, 'ses', session, 'echo', options.template_echo);
                end
            end
        end
    else % If there are NOT sessions:
        session = [];
        % Get tasks for current subject. There is always at least one task (e.g. rest) if there is functional BOLD data, so no need to test.
        tasks = bids.query(options.bids_dataset,'tasks','sub', sub);
        % Loop through tasks
        for t = 1:numel(tasks)
            task = tasks{t};
            % Get runs for current task
            runs = bids.query(options.bids_dataset,'runs','sub', sub, 'task', task);
            % If there are runs, loop through them
            if ~isempty(runs)
                for r = 1:numel(runs)
                    current_run = runs{r};
                    % STEP 1 -- Basic functional preprocessing: fmrwhy_bids_preprocFunc.m
                    fmrwhy_bids_preprocFunc(options.bids_dir, sub, task, options, 'run', current_run);
                    % STEP 2 -- Quality control pipeline: fmrwhy_bids_qcRun.m
                    fmrwhy_bids_qcRun(options.bids_dir, sub, task, options, 'run', current_run, 'echo', options.template_echo);
                end
            else % If there are NOT runs
                % STEP 1 -- Basic functional preprocessing: fmrwhy_bids_preprocFunc.m
                fmrwhy_bids_preprocFunc(options.bids_dir, sub, task, options);
                % STEP 2 -- Quality control pipeline: fmrwhy_bids_qcRun.m
                fmrwhy_bids_qcRun(options.bids_dir, sub, task, options, 'echo', options.template_echo);
            end 
        end
    end
 
    % -------
    % STEP 5 -- QC report: fmrwhy_workflow_qcSubReport.m
    % -------
    fmrwhy_workflow_qcSubReport(sub, options);

end