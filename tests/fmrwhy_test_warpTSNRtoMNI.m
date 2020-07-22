% -------
% STEP 1 -- Load defaults, filenames and parameters
% -------
% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
bids_dir = '/Volumes/Stephan_WD/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

options.dash_me_dir = '/Users/jheunis/Documents/Websites/rt-me-fmri-dash/bids/derivatives/fmrwhy-multiecho';
options.dash_deriv_dir = '/Users/jheunis/Documents/Websites/rt-me-fmri-dash/bids/derivatives';

% Set subject, sessions
subs = {'002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
%subs = {'001'};
tasks = {'rest', 'motor', 'emotion'};
runs = {'1', '2'};
echo = '2';
ses = '';



for s = 1:numel(subs)
    sub = subs{s};

    % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

    % More dirs
    options.me_dir = fullfile(options.deriv_dir, 'fmrwhy-multiecho');
    options.sub_dir_me = fullfile(options.me_dir, ['sub-' sub]);
    options.func_dir_me = fullfile(options.sub_dir_me, 'func');
    dash_sub_dir = fullfile(options.dash_me_dir, ['sub-' sub]);
    if ~exist(dash_sub_dir, 'dir')
        mkdir(dash_sub_dir)
    end

    % Initialise vars
    toTransform_fns = {};
    saveAs_fns = {};
    transformation_fn = [];
    template_fn = [];
    count = 0;

    % Loop through sessions, tasks, runs, etc
    for t = 1:numel(tasks)

        task = tasks{t};

        for r = 1:numel(runs)
            run = runs{r};

            if strcmp(task, 'rest') == 1 && strcmp(run, '1') == 1
                disp('------------')
                disp(['... Skipping Task: ' task ';  Run: ' run ' ...'])
                disp('------------')
                continue;
            end

            disp('------------')
            disp('------------')
            disp(['Task: ' task ';  Run: ' run])
            disp('------------')
            disp('------------')

            % Filenames
            options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

            % Template functional volume
            template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);

            % Transformation from individual functional template to mni space
            transformation_fn = options.indiv_to_mni_fn;

            % Grab tSNR files for each timeseries
            rafunctional_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_tsnr.nii']);
            combined_t2s_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEt2star_tsnr.nii']);
            combined_tsnr_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEtsnr_tsnr.nii']);
            combined_te_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEte_tsnr.nii']);
            tsnr_fns = {rafunctional_fn, combined_t2s_fn, combined_tsnr_fn, combined_te_fn};
            toTransform_fns = [toTransform_fns, tsnr_fns];

            % Create saveAs_fns
            save_rafunctional_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_echo-2_space-MNI152_desc-rapreproc_tsnr.nii']);
            save_combined_t2s_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_space-MNI152_desc-combinedMEt2star_tsnr.nii']);
            save_combined_tsnr_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_space-MNI152_desc-combinedMEtsnr_tsnr.nii']);
            save_combined_te_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_space-MNI152_desc-combinedMEte_tsnr.nii']);
            new_saveAs_fns = {save_rafunctional_fn, save_combined_t2s_fn, save_combined_tsnr_fn, save_combined_te_fn};
            saveAs_fns = [saveAs_fns, new_saveAs_fns];

        end
    end
    rafunctional1_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-rest_run-1_echo-1_desc-rapreproc_tsnr.nii']);
    rafunctional2_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-rest_run-1_echo-2_desc-rapreproc_tsnr.nii']);
    rafunctional3_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-rest_run-1_echo-3_desc-rapreproc_tsnr.nii']);
    run1_tsnr_fns = {rafunctional1_fn, rafunctional2_fn, rafunctional3_fn};
    toTransform_fns = [toTransform_fns, run1_tsnr_fns];

    save_rafunctional1_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-rest_run-1_echo-1_space-MNI152_desc-rapreproc_tsnr.nii']);
    save_rafunctional2_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-rest_run-1_echo-2_space-MNI152_desc-rapreproc_tsnr.nii']);
    save_rafunctional3_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-rest_run-1_echo-3_space-MNI152_desc-rapreproc_tsnr.nii']);
    run1_saveAs_fns = {save_rafunctional1_fn, save_rafunctional2_fn, save_rafunctional3_fn};
    saveAs_fns = [saveAs_fns, run1_saveAs_fns];

    fmrwhy_batch_normaliseWrite(toTransform_fns, transformation_fn, template_fn, saveAs_fns)
end