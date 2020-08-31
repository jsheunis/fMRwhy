% A custom workflow that does ...

% Code steps:
% 1. Define template/default variables, directories and filenames
% 2. Create functional template, if it does not exist
% 3. Estimate 3D volume realignment parameters from raw template echo timeseries (given supplied template volume)
% 4. Run slice time correction for each echo timeseries
% 5. Realign each echo timeseries by applying rigid body transormation estimated from template echo realignment parameters
% 6. Calculate tSNR per echo, using the slice time corrected and realigned functional timeseries as inputs
% 7. Estimate T2star and S0 maps from minimally preprocessed multi-echo data

%--------------------------------------------------------------------------


% -------
% STEP 0.1 -- Load defaults, filenames and parameters
% -------

% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
%bids_dir = '/Volumes/Stephan_WD/NEUFEPME_data_BIDS';
%bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/';
bids_dir = '/Users/jheunis/Desktop/NEUFEPME_data_BIDS/';
%bids_dir = '/Volumes/TSM/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);
options.me_dir = fullfile(options.deriv_dir, 'fmrwhy-multiecho');
if ~exist(options.me_dir, 'dir')
    mkdir(options.me_dir);
end

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Loop through subjects, sessions, tasks, runs, etc
subs = {'001', '002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
subs = {'011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
ses = '';


for s = 1:numel(subs)
    sub = subs{s};
    % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

    % Update workflow params with ME directoires
    options.sub_dir_me = fullfile(options.me_dir, ['sub-' sub]);
    if ~exist(options.sub_dir_me, 'dir')
        mkdir(options.sub_dir_me);
    end
    options.func_dir_me = fullfile(options.sub_dir_me, 'func');
    if ~exist(options.func_dir_me, 'dir')
        mkdir(options.func_dir_me);
    end

    % --------
    % --------
    % STEP 1: Create functional template
    % --------
    % --------
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

    % --------
    % --------
    % STEP 2: For all tasks and runs, complete minimal preprocessing for multi-echo combination
    % --------
    % --------
    tasks = {'rest', 'motor', 'emotion'};
    runs = {'1', '2'};
    for t = 1:numel(tasks)
        task = tasks{t};
        for r = 1:numel(runs)
            run = runs{r};
            fmrwhy_preproc_minFunc(bids_dir, sub, ses, task, run, options)
        end
    end

    % --------
    % --------
    % STEP 3: Run template process on specific task and run predefined the multi-echo template
    % --------
    % --------
    % Template info
    task = 'rest';
    run = '1';

    % --------
    % STEP 3.1: Calculate tSNR per echo, using the slice time corrected and realigned functional timeseries!!!
    % --------
    % load masks (assumes this step has been completed before)
    masks = fmrwhy_util_loadMasks(bids_dir, sub);
    mask_fn = masks.brain_mask_fn;

    for e = 1:options.Ne
        disp('---')
        disp(['Echo ' num2str(e)])
        disp('---')
        % Update workflow params with subject functional derivative filenames
        options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, num2str(e), options);
        tsnr_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_echo-' num2str(e) '_desc-rapreproc_tsnr.nii']);

        if ~exist(tsnr_fn, 'file')
            disp('---')
            disp(['Temporal SNR file does not exist yet. Calculating now...']);
            tsnr_output = fmrwhy_util_calculateTSNR(options.rafunctional_fn, mask_fn, tsnr_fn, options.template_fn)
            disp('Complete!')
            disp('---')
        else
            disp(['Temporal SNR already calculated. Loading now...'])
            tsnr_output = struct;
            nii = nii_tool('load', tsnr_fn);
            tsnr_output.data_3D_tsnr = nii.img;
            disp('---')
        end
    end

    % -------
    % STEP 3.2: Calculate/estimate T2star and S0 maps
    % -------
    t2star_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_t2star.nii']);
    s0_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_s0.nii']);
    if ~exist(t2star_fn, 'file')
        disp('---')
        disp(['T2star map file does not exist yet. Calculating now...']);
        me_fns = {};
        for e = 1:options.Ne
            % Update workflow params with subject functional derivative filenames
            options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, num2str(e), options);
            me_fns{e} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' num2str(e) '_desc-rapreproc_bold.nii']);
        end
        mask_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-brain_mask.nii']);
        MEparams = fmrwhy_util_estimateMEparams(me_fns, options.TE, mask_fn, options.template_fn, t2star_fn, s0_fn);
        disp('Complete!')
        disp('---')
    else
        disp(['T2star map already calculated. Loading now...'])
        MEparams = struct;
        nii = nii_tool('load', t2star_fn);
        MEparams.T2star_3D_thresholded = nii.img;
        nii = nii_tool('load', s0_fn);
        MEparams.S0_3D_thresholded = nii.img;
        disp('---')
    end

    % --------
    % --------
    % STEP 4: Run combination process on all tasks and runs except for the template
    % --------
    % --------
    tasks = {'rest', 'motor', 'emotion'};
    runs = {'1', '2'};

    for t = 1:numel(tasks)
        task = tasks{t};
        for r = 1:numel(runs)
            run = runs{r};

            % Skip template task and run
            if strcmp(task, 'rest') == 1 && strcmp(run, '1') == 1
                disp('------------')
                disp(['Skipping Template files: Task = ' task ';  Run = ' run])
                disp('------------')
                continue;
            end

            % -------
            % STEP 4.1.1: Prepare template data and run multi-echo combination functions
            % -------
            disp('------------')
            disp(['Task: ' task ';  Run: ' run])
            disp('------------')
            % Filenames of combined echo timeseries
            combined_t2s_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEt2star_bold.nii']);
            combined_tsnr_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEtsnr_bold.nii']);
            combined_te_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEte_bold.nii']);
            combined_fns = {combined_t2s_fn, combined_tsnr_fn, combined_te_fn};
            run_combine = 0;
            % If at least one of the files do not exist, run the combination routine
            for x = numel(combined_fns)
                if ~exist(combined_fns{x}, 'file')
                    run_combine = 1;
                else
                    disp(['Combined timeseries exists: ' combined_fns{x}])
                end
            end
            if run_combine
                % TODO: this whole section uses spm_vol instead of dicm2nii, need to fix!!!!
                % Grab and construct parameters (data and weights) for multi-echo combination
                TE = options.TE;
                template_spm = spm_vol(options.template_fn);
                template_dim = template_spm.dim;
                Nt = options.Nscans;
                sz = [template_dim Nt numel(TE)];
                func_data = zeros(sz);
                % Concatenating functional data
                for e = 1:numel(TE)
                    echo = num2str(e);
                    rafunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_bold.nii'])
                    func_data(:,:,:,:,e) = spm_read_vols(spm_vol(rafunctional_fn));
                end
                % Loading weight images
                t2star_img = spm_read_vols(spm_vol(t2star_fn));
                tsnr_data = zeros([template_dim numel(TE)]);
                for e = 1:numel(TE)
                    echo = num2str(e);
                    tsnr_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_echo-' echo '_desc-rapreproc_tsnr.nii']);
                    tsnr_data(:,:,:,e) = spm_read_vols(spm_vol(tsnr_fn));
                end
                % Mult-echo combination
                disp('Combining multiple echoes with different combination methods')
                combined_dataAll_t2s = fmrwhy_me_combineEchoes(func_data, TE, 0, 1, t2star_img);
                combined_dataAll_tsnr = fmrwhy_me_combineEchoes(func_data, TE, 0, 2, tsnr_data);
                combined_dataAll_TE = fmrwhy_me_combineEchoes(func_data, TE, 0, 3, TE);

                % Save nifti images
                % TODO: replace this with dicm2nii equivalent (ask on github repo)
                rafunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_bold.nii']);
                new_spm_t2s = spm_vol(rafunctional_fn);
                new_spm_tsnr = new_spm_t2s;
                new_spm_te = new_spm_t2s;

                for i = 1:numel(new_spm_t2s)
                    new_spm_t2s(i).fname = combined_t2s_fn;
                    new_spm_t2s(i).private.dat.fname = combined_t2s_fn;
                    spm_write_vol(new_spm_t2s(i), combined_dataAll_t2s(:,:,:,i));

                    new_spm_tsnr(i).fname = combined_tsnr_fn;
                    new_spm_tsnr(i).private.dat.fname = combined_tsnr_fn;
                    spm_write_vol(new_spm_tsnr(i), combined_dataAll_tsnr(:,:,:,i));

                    new_spm_te(i).fname = combined_te_fn;
                    new_spm_te(i).private.dat.fname = combined_te_fn;
                    spm_write_vol(new_spm_te(i), combined_dataAll_TE(:,:,:,i));
                end
            end

            % -------
            % STEP 4.1.2: Prepare template data and run FIT multi-echo combination
            % -------
            combined_t2sFIT_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEt2starFIT_bold.nii']);
            t2sFIT_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-t2starFIT_bold.nii']);
            s0FIT_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-s0FIT_bold.nii']);
            fit_fns = {combined_t2sFIT_fn, t2sFIT_fn, s0FIT_fn};
            run_combine = 0;
            % If at least one of the files do not exist, run the combination routine
            for x = numel(fit_fns)
                if ~exist(fit_fns{x}, 'file')
                    run_combine = 1;
                else
                    disp(['Timeseries exists: ' fit_fns{x}])
                end
            end
            if run_combine
                TE = options.TE;
                template_spm = spm_vol(options.template_fn);
                template_dim = template_spm.dim;
                Nt = options.Nscans;
                sz = [template_dim Nt numel(TE)];
                func_data = zeros(sz);
                % Concatenating functional data
                for e = 1:numel(TE)
                    echo = num2str(e);
                    rafunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_bold.nii']);
                    func_data(:,:,:,:,e) = spm_read_vols(spm_vol(rafunctional_fn));
                end

                % TODO: replace this with dicm2nii equivalent (ask on github repo)
                rafunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_bold.nii']);
                new_spm_combined_t2sFIT = spm_vol(rafunctional_fn);
                new_spm_t2sFIT = new_spm_combined_t2sFIT;
                new_spm_s0FIT = new_spm_combined_t2sFIT;
                masksSPM = fmrwhy_util_loadMasksSPM(bids_dir, sub);
                I_mask = masksSPM.brain_mask_I;

                for i = 1:Nt

                    func_data_pv = {};
                    for e = 1:numel(TE)
                        func_data_pv{e} = squeeze(func_data(:,:,:,i,e));
                    end

                    me_params = fmrwhy_realtime_estimateMEparams(func_data_pv, TE, I_mask);

                    func_pv = zeros([template_dim numel(TE)]);
                    for e = 1:numel(TE)
                        func_pv(:,:,:,e) = func_data_pv{e};
                    end

                    combined_t2sFIT_3D = fmrwhy_me_combineEchoes(func_pv, TE, 0, 1, me_params.T2star_3D_thresholded);

                    new_spm_combined_t2sFIT(i).fname = combined_t2sFIT_fn;
                    new_spm_combined_t2sFIT(i).private.dat.fname = combined_t2sFIT_fn;
                    spm_write_vol(new_spm_combined_t2sFIT(i), combined_t2sFIT_3D);

                    new_spm_t2sFIT(i).fname = t2sFIT_fn;
                    new_spm_t2sFIT(i).private.dat.fname = t2sFIT_fn;
                    new_spm_t2sFIT(i).pinfo(1) = 1;
                    spm_write_vol(new_spm_t2sFIT(i), me_params.T2star_3D_thresholded);

                    new_spm_s0FIT(i).fname = s0FIT_fn;
                    new_spm_s0FIT(i).private.dat.fname = s0FIT_fn;
                    spm_write_vol(new_spm_s0FIT(i), me_params.S0_3D_thresholded);

                end
            end



            % -------
            % STEP 4.2: Calculate tSNR for each combined timeseries
            % -------
            rafunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_bold.nii']);
            main_fns = {rafunctional_fn, combined_t2s_fn, combined_tsnr_fn, combined_te_fn};
            main_fns = [main_fns, fit_fns];
            tsnr_fns = {};
            tsnr_output = {};
            for i = 1:numel(main_fns)
                if i == 1
                    tsnr_fns{i} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_tsnr.nii'])
                else
                    tsnr_fns{i} = strrep(main_fns{i}, 'bold', 'tsnr');
                end

                if ~exist(tsnr_fns{i}, 'file')
                    tsnr_output{i} = fmrwhy_util_calculateTSNR(main_fns{i}, mask_fn, tsnr_fns{i}, template_fn);
                else
                    disp(['tSNR already computed for timeseries: ' tsnr_fns{i}])
                end
            end

            % -------
            % STEP 4.3: Smooth each combined timeseries, for later analysis purposes
            % -------
            smooth_fns = {};
            for i = 1:numel(main_fns)
                smooth_fns{i} = strrep(main_fns{i}, 'desc-', 'desc-s')
                if (~exist(smooth_fns{i}, 'file')) && (i<7)
                    fmrwhy_batch_smooth(main_fns{i}, smooth_fns{i}, options.fwhm);
                else
                    disp(['Spatial smoothing already completed for timeseries: ' smooth_fns{i}])
                end
            end

        end
    end
end
