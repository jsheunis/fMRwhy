% -------
% STEP 1 -- Load defaults, filenames and parameters
% -------
% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Set subject, sessions
sub = '001';
ses = '';

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% load mask
masks = fmrwhy_util_loadMasks(bids_dir, sub);
mask_fn = masks.brain_mask_fn;

% Loop through sessions, tasks, runs, etc
tasks = {'rest', 'motor', 'emotion'};
runs = {'1', '2'};
echo = '2';

for t = 1:numel(tasks)

    task = tasks{t};

    for r = 1:numel(runs)
        run = runs{r};

        if strcmp(task, 'rest') == 1 && strcmp(run, '1') == 1
            disp('------------');
            disp(['... Skipping Task: ' task ';  Run: ' run ' ...']);
            disp('------------');
            continue
        end

        disp('------------');
        disp('------------');
        disp(['Task: ' task ';  Run: ' run]);
        disp('------------');
        disp('------------');

        % Filenames
        options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

        % Calculate tSNR for each timeseries
        rafunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_bold.nii']);
        combined_t2s_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEt2star_bold.nii']);
        combined_tsnr_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEtsnr_bold.nii']);
        combined_te_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEte_bold.nii']);
        main_fns = {rafunctional_fn, combined_t2s_fn, combined_tsnr_fn, combined_te_fn};
        tsnr_fns = {};
        tsnr_output = {};
        for i = 1:numel(main_fns)
            tsnr_fns{i} = strrep(main_fns{i}, 'bold', 'tsnr');
            tsnr_output{i} = fmrwhy_util_calculateTSNR(main_fns{i}, mask_fn, tsnr_fns{i}, template_fn);
        end
    end
end
