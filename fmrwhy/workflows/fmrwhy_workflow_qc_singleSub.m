function fmrwhy_workflow_qc_singleSub(sub, sessions, tasks, runs, settings_fn, spm_dir)

    % -----------
    % DESCRIPTION
    % -----------
    % An SPM12- and fMRwhy-based preprocessing and QC workflow for a single subject.
    % It contains the following  steps in order:

    % 1. Slice timing correction on all functional runs
    % 2. Realign (estimate + write) all volumes of all functional runs to the mean functional image (across runs).
    % This is a three-step process:
    %   - First realign first volumes of all runs to the first volume of the first run
    %   - Then realign all volumes per run to the first volume of the applicable run
    %   - Then calculate the mean image across runs and realign all volumes of all runs to this mean image
    % 3. Coregister (estimate) the anatomical image to the mean functional volume. Only estimate transform and write to image header,
    %    do not reslice yet. The anatomical image is now in functional space, but not yet functional resolution.
    % 4. Normalise (estimate + write) the coregistered anatomical image to MNI space.
    %    Apply the same transform to all realigned functional runs. All images are now in standard MNI space.
    %    Reslice/resample all to specified resolution (2x2x2mm).
    % 5. Spatial smoothing of all functional runs

    % ---------------------
    % SOFTWARE DEPENDENCIES
    % ---------------------
    % 1. fMRwhy (lcid branch): https://github.com/jsheunis/fMRwhy/tree/lcid
    % 2. SPM12:  https://github.com/spm/spm12/releases/tag/r7771
    % 3. bids-matlab:  https://github.com/bids-standard/bids-matlab
    % 4. dicm2nii:  https://github.com/jsheunis/dicm2nii/releases/tag/v0.2
    % 5. TAPAS PhysIO:  https://github.com/translationalneuromodeling/tapas/releases/tag/v4.0.0
    % 6. Raincoud Plots:  https://github.com/RainCloudPlots/RainCloudPlots/releases/tag/v1.1

    % -------------------
    % SETUP (DO NOT EDIT)
    % -------------------
    % SETUP STEP A -- Set fMRwhy defaults, check dependencies, set Matlab path, etc
    options = fmrwhy_defaults();
    options = fmrwhy_util_checkDependencies(options);
    % SETUP STEP B -- Load settings, directories, filenames and parameters
    % Run settings file ==> populates study/data-specific fields in the options structure, including BIDS variables
    run(settings_fn);
    % Setup fmrwhy derivative directories on workflow level
    options = fmrwhy_bids_setupQcDerivDirs(options.bids_dir, options);
    % Validate settings
    options = fmrwhy_settings_validate(options);

    % -------------------------------------------
    % SUBJECTS, SESSIONS, TASKS, RUNS (DO NOT EDIT THIS, unless you specifically want to override defaults)
    % -------------------------------------------
    % Create a cell array with the required sessions to run
    % sessions = {'1', '2', '3'};
    % Create a cell array with the required tasks
    % tasks = {'rest', 'motor'};
    % Create a cell array with the required runs
    % runs = {'1', '2', '3'};

    % NOTE: for now, we're assuming that all subjects have all sessions, all sessions have all tasks, and all tasks has all runs.
    % During the workflow, we loop through all possible combinations of these,
    % and we do a check to see if the corresponding files actually exist before starting the run-specific pipeline.
    % If it doesn't exist, we write a notification to the command line and then continue with the next iteration.
    % There are smarter ways of building this logic, but this is the easiest for now. TODO.

    % ----------------------
    % PREPROCESSING PIPELINE
    % ----------------------

    % The following is executed per session, and per task.
    % For runs, the preprocessing steps are either executed on the task level (i.e. all runs in one step),
    % or on the individual run level, depending on the specific preprocessing step.
    % TODO: this currently assumes that there are sessions, need to fully address the possibility of no sessions.

    % Setup fmrwhy derivatives directories on subject level (this copies data from the main bids_dir)
    options = fmrwhy_bids_setupQcSubDirs(options.bids_dir, sub, options);

    for ss = 1:numel(sessions)
        ses = sessions{ss};

        % Skip session loop if current session does not exist
        ses_dir = fullfile(options.preproc_dir, ['sub-' sub], ['ses-' ses]);
        if ~exist(ses_dir, 'dir')
            disp('---');
            disp(['No session ses-', ses, ' for sub-', sub]);
            disp('---');
            continue
        end

        % Skip session loop if there is no functional data for the current session
        func_dir = fullfile(options.preproc_dir, ['sub-' sub], ['ses-' ses], 'func');
        if ~exist(func_dir, 'dir')
            disp('---');
            disp(['No functional data for sub-', sub, '_ses-', ses]);
            disp('---');
            continue
        end

        % Set up anatomical data for coregistration
        % If a template session is specified and it is not equal to the current session,
        % copy the template data to current session
        options = fmrwhy_bids_getAnatDerivs(options.bids_dir, sub, options, 'ses', ses);
        if isempty(options.anat_template_session) && (ses ~= options.anat_template_session)
            [fn, fp] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'ses', options.anat_template_session, 'ext', '_T1w.nii');
            anatomical_fn = fullfile(options.preproc_dir, fp, fn);
            [fn_new, fp_new] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'ses', ses, 'ext', '_T1w.nii');
            new_anatomical_fn = fullfile(options.preproc_dir, fp_new, fn_new);
            copyfile(anatomical_fn, new_anatomical_fn);
            options.anat_template_session = '';
        end

        % Loop through tasks
        for t = 1:numel(tasks)
            task = tasks{t};

            % Check if files with current task are available; if not, skip iteration
            fnames = dir([func_dir filesep '*' task '*']);
            if isempty(fnames)
                disp('---');
                disp(['No task files for sub-', sub, '_ses-', ses, '_task-', task]);
                disp('---');
                continue
            end

            disp('---------------------');
            disp(['START PREPROCESSING: sub-', sub, '_ses-', ses, '_task-', task]);
            disp('---------------------');

            % -------
            % STEP 0 -- Calculate realignment parameters per run (important: before slice timing correction)
            % -------
            for r = 1:numel(runs)
                rn = runs{r};
                % Update functional derivate filenames
                options = fmrwhy_bids_getFuncDerivs(options.bids_dir, sub, task, options, 'ses', ses, 'task', task, 'run', rn);
                if exist(options.functional_fn, 'file')
                    disp('---');
                    disp(['STEP 0) Calculate realignment parameters: sub-' sub '_ses-' ses '_task-' task '_run-' rn]);
                    disp('---');
                    % Check if realignment has already been done by seeing if the tsv file with head movement parameters exist
                    [d, f, e] = fileparts(options.motion_fn);
                    if ~exist(options.motion_fn, 'file')
                        % If it does not exist estimate MPs
                        realign_measures = fmrwhy_batch_realignEst(options.functional_fn, 0); % 0 ==> use first volume of each run as run's template for estimation
                        temp_txt_fn = fullfile(d, [f '.txt']);
                        col_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
                        data = load(realign_measures.mp_fn);
                        data_table = array2table(data, 'VariableNames', col_names);
                        writetable(data_table, temp_txt_fn, 'Delimiter', '\t');
                        [status, msg, msgID] = movefile(temp_txt_fn, options.motion_fn);
                    else
                        disp(['3D realignment parameters already estimated: ' options.motion_fn]);
                        disp('---');
                    end
                else
                    disp(['Functional file not available: ' options.functional_fn]);
                end
            end

            % -------
            % STEP 1 -- Slice timing correction
            % -------
            for r = 1:numel(runs)
                rn = runs{r};
                options = fmrwhy_bids_getFuncDerivs(options.bids_dir, sub, task, options, 'ses', ses, 'task', task, 'run', rn);
                if exist(options.functional_fn, 'file')
                    disp('---');
                    disp(['STEP 1) Slice Timing Correction: sub-' sub '_ses-' ses '_task-' task '_run-' rn]);
                    disp('---');
                    if exist(options.afunctional_fn, 'file')
                        disp(['slice timing correction already completed...']);
                        disp('---');
                    else
                        fmrwhy_batch_sliceTiming(options.functional_fn, options.afunctional_fn, options);
                    end
                else
                    disp(['Functional file not available: ' options.functional_fn]);
                end
            end

            % -------
            % STEP 2 -- Realignment
            % -------
            switch options.realignment_type
                case 'per_task'
                    % All runs of a task are included (in order) in the realignment procedure; a two step procedure is applied.
                    afunctional_fns = {};
                    saveAs_fns = {};
                    for r = 1:numel(runs)
                        rn = runs{r};
                        options = fmrwhy_bids_getFuncDerivs(options.bids_dir, sub, task, options, 'ses', ses, 'task', task, 'run', rn);
                        if exist(options.functional_fn, 'file')
                            afunctional_fns = [afunctional_fns {options.afunctional_fn}];
                            saveAs_fns = [saveAs_fns {options.rafunctional_fn}];
                        else
                            disp(['Functional file not available: ' options.functional_fn]);
                        end
                    end
                    disp('---');
                    disp(['STEP 2) 3D Volume Realignment: sub-' sub '_ses-' ses '_task-' task '_run-all']);
                    disp('---');
                    if exist(saveAs_fns{1}, 'file')
                        disp(['3D volume realignment already completed...']);
                        disp('---');
                    else
                        fmrwhy_batch_realignEstResl(afunctional_fns, 0, saveAs_fns, 'fwhm', 6, 'rtm', 1, 'which', [2 1]);
                        % TODO, this does not handle the mean image output renaming in any standard way yet. Standardise.
                    end
                case 'per_run'
                    disp('per_run');
                case 'to_template'
                    disp('to_template');
                otherwise
                    warning('Unexpected realignmet type. realignment skipped.');
            end

            % -------
            % STEP 3 -- Structural-functional preprocessing
            % Including:
            % - Anatomical to functional space coregistration - SPM12 coregister estimate
            % - Segment coregistered anatomical image into tissue components - SPM12 unified segmentation
            %     - (this saves inverse transform from subject functional to MNI space)
            % - Reslice all to functional space grid (SPM reslice)
            % - Create tissue compartment and whole brain masks
            % -------
            disp('---');
            disp(['STEP 3) Structural-functional preprocessing: sub-' sub '_ses-' ses]);
            disp('---');
            switch options.coreg_type
                case 'per_task'
                    % Get anatomical derivative filenames
                    options = fmrwhy_bids_getAnatDerivs(options.bids_dir, sub, options, 'ses', ses, 'task', task);
                    % Add template filename (mean image from realignment)
                    [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', sub, 'ses', ses, 'task', task, 'run', runs{1}, 'desc', 'apreproc', 'ext', '_bold.nii');
                    options.template_fn = fullfile(options.preproc_dir, filepath, ['meantemp_' filename]);
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
                        disp('Running complete structural-funcional preprocessing pipeline');
                        fmrwhy_preproc_structFunc(options);
                        disp('Complete!');
                        disp('---');
                    else
                        disp('Structural-funcional preprocessing already completed.');
                        disp('---');
                    end
                case 'per_run'
                    disp('per_run');
                case 'to_template'
                    disp('to_template');
                otherwise
                    warning('Unexpected coregistration type. Coregistration skipped.');
            end

            % -------
            % STEP 4 -- Normalise (estimate+write) all
            % -------
            toTransform_fns = {};
            saveAs_fns = {};
            for r = 1:numel(runs)
                rn = runs{r};
                options = fmrwhy_bids_getFuncDerivs(options.bids_dir, sub, task, options, 'ses', ses, 'task', task, 'run', rn);
                if exist(options.rafunctional_fn, 'file')
                    toTransform_fns = [toTransform_fns {options.rafunctional_fn}];
                    % TODO check/get name of warped file in MNI space
                    [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', sub, 'ses', ses, 'task', task, 'run', rn, 'desc', 'wrapreproc', 'ext', '_bold.nii');
                    % Probably with fmrwhy_bids_constructFilename
                    saveAs_fns = [saveAs_fns {fullfile(options.preproc_dir, filepath, filename)}];
                else
                    disp(['WARNING -- Realigned and slice-timing-corrected functional file not found: ' options.rafunctional_fn]);
                end
            end
            disp('---');
            disp(['STEP 4) Normalise to MNI space: sub-' sub '_ses-' ses '_task-' task '_run-all']);
            disp('---');
            reference_fn = options.coregest_anatomical_fn;
            tpm_fn = fullfile(spm_dir, 'tpm/TPM.nii');
            if exist(saveAs_fns{1}, 'file')
                disp(['normalisation already completed...']);
                disp('---');
            else
                fmrwhy_batch_normaliseEstWrite(reference_fn, toTransform_fns, tpm_fn, saveAs_fns);
            end

            % -------
            % STEP 5 -- Spatial smoothing
            % -------
            disp('---');
            disp('STEP 5) Spatial Smoothing');
            disp('---');
            for r = 1:numel(runs)
                rn = runs{r};
                [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', sub, 'ses', ses, 'task', task, 'run', rn, 'desc', 'wrapreproc', 'ext', '_bold.nii');
                wrafunctional_fn = fullfile(options.preproc_dir, filepath, filename);
                [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', sub, 'ses', ses, 'task', task, 'run', rn, 'desc', 'swrapreproc', 'ext', '_bold.nii');
                swrafunctional_fn = fullfile(options.preproc_dir, filepath, filename);
                if exist(wrafunctional_fn, 'file')
                    if exist(swrafunctional_fn, 'file')
                        disp(['spatial smoothing already completed...']);
                        disp('---');
                    else
                        disp(['Spatial Smoothing: sub-' sub '_ses-' ses '_task-' task '_run-' rn]);
                        fmrwhy_batch_smooth(wrafunctional_fn, swrafunctional_fn, options.fwhm);
                    end
                else
                    disp(['WARNING -- Warped, realigned and slice-timing-corrected functional file not found: ' wrafunctional_fn]);
                end
            end

            % -------
            % QUALITY CONTROL STEPS
            % -------
            disp('---');
            disp('STEP 6) Quality control processing');
            disp('---');
            % -------
            % A - smooth raw functional data for the plot;
            disp('STEP 6-A) Smooth raw functional data for the plot');
            disp('---');
            for r = 1:numel(runs)
                rn = runs{r};
                options = fmrwhy_bids_getFuncDerivs(options.bids_dir, sub, task, options, 'ses', ses, 'task', task, 'run', rn);
                % Smooth raw timeseries data
                if exist(options.functional_fn, 'file')
                    if ~exist(options.sfunctional_fn, 'file')
                        disp(['Performing spatial smoothing on raw timeseries: ' options.current_functional_filename]);
                        fmrwhy_batch_smooth(options.functional_fn, options.sfunctional_fn, options.fwhm);
                        disp('Complete!');
                        disp('---');
                    else
                        disp(['Spatial smoothing already completed for raw timeseries: ' options.current_functional_filename]);
                        disp('---');
                    end
                else
                    disp(['WARNING -- Raw functional file not found for run: ' options.functional_fn]);
                end
            end

            % -------
            % B - Generate multiple regressors for GLM analysis and QC
            % Includes:
            % - 3D realignment parameters
            % - framewise displacement
            % - FD censoring
            % - tissue compartment signals
            % -------
            disp('STEP 6-B) Generate multiple regressors for GLM analysis and QC');
            disp('---');
            for r = 1:numel(runs)
                rn = runs{r};
                options = fmrwhy_bids_getFuncDerivs(options.bids_dir, sub, task, options, 'ses', ses, 'task', task, 'run', rn);
                if exist(options.functional_fn, 'file')
                    if ~exist(options.confounds_fn, 'file')
                        disp(['Generating multiple regressors for GLM analysis and QC']);
                        fmrwhy_bids_preprocMultRegr(options.bids_dir, sub, task, options, 'ses', ses, 'task', task, 'run', rn);
                        disp('Complete!');
                        disp('---');
                    else
                        disp(['Multiple regressors already generated: ' options.confounds_fn]);
                        disp('---');
                    end
                else
                    disp(['WARNING -- Functional file not found for run: ' options.functional_fn]);
                end
            end

            % -------
            % C - Calculate QC metrics and generate plots
            % -------
            disp('STEP 6-C) Calculate QC metrics and generate plots');
            disp('---');
            for r = 1:numel(runs)
                rn = runs{r};
                options = fmrwhy_bids_getFuncDerivs(options.bids_dir, sub, task, options, 'ses', ses, 'task', task, 'run', rn);
                if exist(options.functional_fn, 'file')
                    fmrwhy_bids_qcRun(options.bids_dir, sub, task, options, 'ses', ses, 'run', rn);
                else
                    disp(['WARNING -- Functional file not found for run: ' options.functional_fn]);
                end
            end

        end
    end

    % D - Generate subject report
    disp('STEP 7) Generate subject QC report');
    disp('---');
    fmrwhy_workflow_qcSubReport(sub, options);
    disp('Complete!');
    disp('---');

end
