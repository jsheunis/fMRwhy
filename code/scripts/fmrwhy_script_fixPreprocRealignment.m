% -------
% STEP 1 -- Load defaults, filenames and parameters
% -------
% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
%bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/';
bids_dir = '/Users/jheunis/Desktop/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Set subject, sessions
subs = {'016', '017', '018', '019', '020', '021', '022', '023', '024'};
%subs = {'002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
ses = '';


for s = 1:numel(subs)
    tic;
    sub = subs{s};


    % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

    % Loop through sessions, tasks, runs, etc
    tasks = {'rest', 'motor', 'emotion'};
    runs = {'1', '2'};

    for t = 1:numel(tasks)

        task = tasks{t};

        for r = 1:numel(runs)
            run = runs{r};

            disp('------------')
            disp('------------')
            disp(['Task: ' task ';  Run: ' run])
            disp('------------')
            disp('------------')

            % -------
            % STEP 1: First access template timeseries information
            % -------
            options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, options.template_echo, options);

            % Skip slice timing correction of raw timeseries, this should have been implemented correctly previously

            % -------
            % STEP 2: 3D volume realignment
            % -------
            % The first attempt at getting realignment parameters was already successful, so no need to rerun.
            motion_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' options.template_echo '_desc-confounds_motion.tsv']);
            motion_struct = tdfread(motion_fn)
            motion_params = struct2array(motion_struct);

            for e = 1:numel(options.TE)
                disp('---')
                disp(['Echo ' num2str(e)])
                disp('---')
                % Filler text
                stre_txt = ['sub-' sub '_task-' task '_run-' run '_echo-' num2str(e)];
                % Update workflow params with subject functional derivative filenames
                options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, num2str(e), options);
                % Realign raw timeseries data
                disp(['Performing 3D realignment on raw timeseries: ' stre_txt])
                fmrwhy_util_applyTransform(options.functional_fn, motion_params, options.template_fn, options.rfunctional_fn)
                disp('Complete!')
                disp('---')

                % Realign slice time corrected timeseries data
                disp(['Performing 3D realignment on slice time corrected timeseries: ' stre_txt])
                fmrwhy_util_applyTransform(options.afunctional_fn, motion_params, options.template_fn, options.rafunctional_fn)
                disp('Complete!')
                disp('---')
            end
            % -------
            % STEP 3: spatial smoothing
            % -------
            for e = 1:numel(options.TE)
                disp('---')
                disp(['Echo ' num2str(e)])
                disp('---')
                % Filler text
                stre_txt = ['sub-' sub '_task-' task '_run-' run '_echo-' num2str(e)];
                % Update workflow params with subject functional derivative filenames
                options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, num2str(e), options);

                % Smoothing of raw timeseries data already done correctly
                % Smooth realigned timeseries data
                disp(['Performing spatial smoothing on realigned timeseries: ' stre_txt])
                fmrwhy_batch_smooth(options.rfunctional_fn, options.srfunctional_fn, options.fwhm);
                disp('Complete!')
                disp('---')
                % Smooth realigned and slice time corrected timeseries data
                disp(['Performing spatial smoothing on realigned and slice time corrected timeseries: ' stre_txt])
                fmrwhy_batch_smooth(options.rafunctional_fn, options.srafunctional_fn, options.fwhm);
                disp('Complete!')
                disp('---')
            end

            % -------
            % STEP 4: fix tissue compartment signals in GLM regressors file
            % -------
            % Get correct variable references
            echo = options.template_echo;
            options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);
            % First delete the files with incorrect information (i.e. those that used incorrectly realigned data)
            delete(options.tissue_regr_fn)
            delete(options.confounds_fn)
            % Then run the function to re-generate these files
            confounds_tsv = fmrwhy_preproc_generateMultRegr(bids_dir, sub, ses, task, run, echo, options)

            % -------
            % STEP 5: fix QC processing and outputs (including redoing QC image styling)
            % -------
            % Recreating the plot for whole brain and ROI are not necessary, as they use options.sfunctional_fn (i.e. not realigned)
            fmrwhy_qc_run(bids_dir, sub, ses, task, run, options.template_echo, options);

            % -------
            % STEP 6: fix ME processing and outputs
            % -------
%            options.me_dir = fullfile(options.deriv_dir, 'fmrwhy-multiecho');
%            options.sub_dir_me = fullfile(options.me_dir, ['sub-' sub]);
%            options.func_dir_me = fullfile(options.sub_dir_me, 'func');
%            all_files = dir(fullfile(options.func_dir_me, ['sub-' sub '*']));
%            for i=1:numel(all_files)
%                delete(fullfile(options.func_dir_me, all_files(i).name))
%            end
        end
    end
    fmrwhy_neufep_generateSubReport(bids_dir, sub)
    toc;
end