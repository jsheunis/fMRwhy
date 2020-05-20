% A custom workflow that does structFunc and basicFunc preprocessing and QC for a single subject in the NEUFEP study

% Code steps:
% 1. Define template/default variables, directories and filenames
% 2. Specify


%--------------------------------------------------------------------------


% -------
% STEP 0.1 -- Load defaults, filenames and parameters
% -------

% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Loop through subjects, sessions, tasks, runs, etc
subs = {'003', '004', '005'};
%sub = '002';
ses = '';


for s = 1:numel(subs)
    sub = subs{s};
    % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

    % -------
    % STEP 0.2 -- Create functional template
    % -------
    % Create, if it does not exist
    template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
    if ~exist(template_fn, 'file')
        disp(['Template funcional image does not exist yet. Creating now: ' template_fn]);
        functional_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_echo-' options.template_echo '_bold.nii']);
        fmrwhy_util_saveNiftiFrom4D(functional_fn, template_fn, 1)
    else
        disp(['Template functional image exists: ' template_fn]);
    end
    options.template_fn = template_fn;


    % -------
    % STEP 1 -- Structural-functional preprocessing: fmrwhy_preproc_structFunc.m
    % -------
    % Loop through all standard output filenames and see if these files exist
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
        fmrwhy_preproc_structFunc(bids_dir, sub, ses, options.template_task, options.template_run, options.template_echo, options);
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
                    montage_fn = fullfile(options.anat_dir_qc, ['sub-' sub '_space-individual_desc-' options.roi.(options.tasks{i}).desc{j} '_roi_montage.png']);
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
            fmrwhy_preproc_anatLocaliser(bids_dir, sub, options)
            disp('Complete!')
            disp('---')
        else
            disp('Anatomical localiser processing already completed.')
            disp('---')
        end
    end

    % -------
    % PER TASK and RUN
    % -------

    % Loop through sessions, tasks, runs, echoes.
    ses = '';
    tasks = {'rest', 'motor', 'emotion'};
    runs = {'1', '2'};

    for t = 1:numel(tasks)
        task = tasks{t};
        for r = 1:numel(runs)
            run = runs{r};

            % -------
            % STEP 1 -- Basic functional preprocessing: fmrwhy_preproc_basicFunc.m
            % -------
            % NOTE: all outputs for multi-echo are many files, they will be checked individually in fmrwhy_preproc_basicFunc
            fmrwhy_preproc_basicFunc(bids_dir, sub, ses, task, run, options);

            % -------
            % PREPROC STEP 2 -- Quality control pipeline: fmrwhy_qc_run.m
            % -------
            % NOTE: all outputs for multi-echo are many files, they will be checked individually in fmrwhy_qc_run
            fmrwhy_qc_run(bids_dir, sub, ses, task, run, options.template_echo, options);

        end
    end

    % -------
    % STEP 5 -- QC report: fmrwhy_qc_generateSubRunReport.m
    % -------
    fmrwhy_neufep_generateSubReport(bids_dir, sub);

end