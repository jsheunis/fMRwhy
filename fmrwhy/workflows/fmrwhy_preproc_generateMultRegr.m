function confounds_tsv = fmrwhy_preproc_generateMultRegr(bids_dir, sub, ses, task, run, echo, options)
    % For now, this function assumes that all necessary preproc has been done and all necessary files exist
    % Generates multiple regressors for quality control and

    % Info: These are examples of what can be included, they are not yet all included (TODO)
    % trans_x trans_y trans_z rot_x rot_y rot_z framewise_displacement framewise_displacement_censor global_signal white_matter csf
    % dvars std_dvars dvars_censor std_dvars_censor
    % trans_x_power2, trans_y_power2, trans_z_power2, rot_x_power2, rot_y_power2, rot_z_power2, trans_x_derivative1, trans_y_derivative1, trans_z_derivative1, rot_x_derivative1, rot_y_derivative1, rot_z_derivative1, trans_x_derivative1_power2, trans_y_derivative1_power2, trans_z_derivative1_power2, rot_x_derivative1_power2, rot_y_derivative1_power2, rot_z_derivative1_power2
    % retroicor_c1 retroicor_c2 retroicor_c3 retroicor_c4 retroicor_c5 retroicor_c6
    % retroicor_r1 retroicor_r2 retroicor_r3 retroicor_r4 retroicor_r5 retroicor_r6 retroicor_r7 retroicor_r8
    % hrv rvt

    % Info: regressors from PhysIO explained:
    % RETROICOR cardiac regressors [2 x nOrderCardiac] ==> 2 x 3 = 6
    % RETROICOR respiratory regressors [2 x nOrderRespiratory] ==> 2 x 4 = 8
    % RETROICOR cardXResp interaction regressors [4 x nOrderCardiacXRespiratory] ==> 4 x 1 = 4
    % HRV [nDelaysHRV]
    % RVT [nDelaysRVT]

    % Order of regressors (parts in brackets are dependent on user parameters or not yet included):
    % 1. MP: trans_x trans_y trans_z rot_x rot_y rot_z
    % 2: FD: framewise_displacement framewise_displacement_censor02 framewise_displacement_censor05 (framewise_displacement_censor)
    % 3: (DVARS: dvars std_dvars dvars_censor std_dvars_censor)
    % 4. Tissue: global_signal white_matter csf
    % 5. MP volterra: trans_x_power2, trans_y_power2, trans_z_power2, rot_x_power2, rot_y_power2, rot_z_power2, trans_x_derivative1, trans_y_derivative1, trans_z_derivative1, rot_x_derivative1, rot_y_derivative1, rot_z_derivative1, trans_x_derivative1_power2, trans_y_derivative1_power2, trans_z_derivative1_power2, rot_x_derivative1_power2, rot_y_derivative1_power2, rot_z_derivative1_power2
    % 6. RETROICOR cardiac: retroicor_c1 retroicor_c2 retroicor_c3 retroicor_c4 retroicor_c5 retroicor_c6
    % 7. RETROICOR respiratory: retroicor_r1 retroicor_r2 retroicor_r3 retroicor_r4 retroicor_r5 retroicor_r6 retroicor_r7 retroicor_r8
    % 8. RETROICOR interaction: retroicor_cxr1 retroicor_cxr2 retroicor_cxr3 retroicor_cxr4
    % 9. hrv rvt

    % ---

    % Setup fmrwhy BIDS-derivatuve directories on workflow level
    options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

    % Grab parameters from workflow settings file
    options = fmrwhy_settings_preprocQC(bids_dir, options);

    % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

    % Update workflow params with subject functional derivative filenames
    options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

    % Filler text
    stre_txt = ['sub-' sub '_task-' task '_run-' run '_echo-' echo];
    str_txt = ['sub-' sub '_task-' task '_run-' run];

    % Default steps to include/exclude
    include_volterra = options.confounds.include_volterra;
    include_fd = options.confounds.include_fd;
    include_tissue = options.confounds.include_tissue;
    include_physio = options.confounds.include_physio;

    % -------
    % STEP 1: Generate realignment, framewise displacement and censoring regressors
    % -------
    % First realignment parameters
    % Check if this has already been done by seeing if the tsv file with head movement parameters exist
    [d, f, e] = fileparts(options.motion_fn);
    if ~exist(options.motion_fn, 'file')
        % If it does not exist, estimate MPs
        disp(['Estimating 3D realignment parameters for: ' stre_txt]);
        realign_measures = fmrwhy_batch_realignEst(options.functional_fn, options.template_fn);
        temp_txt_fn = fullfile(d, [f '.txt']);
        col_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
        MP = realign_measures.MP;
        data_table = array2table(MP, 'VariableNames', col_names);
        writetable(data_table, temp_txt_fn, 'Delimiter', '\t');
        [status, msg, msgID] = movefile(temp_txt_fn, options.motion_fn);
        disp('Complete!');
        disp('---');
    else
        disp(['3D realignment parameters already estimated: ' options.motion_fn]);
        disp('---');
        MP = struct2array(tdfread(options.motion_fn));
    end

    % Then do volterra expansion of motion parameters
    % Check if this has already been done by seeing if the tsv file with head movement parameters exist
    if include_volterra
        MP_volterra = fmrwhy_util_createVolterraExpansion(MP, 0);
        col_names_volterra = {'trans_x_power2', 'trans_y_power2', 'trans_z_power2', 'rot_x_power2', 'rot_y_power2', 'rot_z_power2', 'trans_x_derivative1', 'trans_y_derivative1', 'trans_z_derivative1', 'rot_x_derivative1', 'rot_y_derivative1', 'rot_z_derivative1', 'trans_x_derivative1_power2', 'trans_y_derivative1_power2', 'trans_z_derivative1_power2', 'rot_x_derivative1_power2', 'rot_y_derivative1_power2', 'rot_z_derivative1_power2'};
    end

    % Then framewise d isplacement and censoring
    if include_fd
        [d, f, e] = fileparts(options.framewise_displacement_fn);
        if ~exist(options.framewise_displacement_fn, 'file')
            % If it does not exist, estimate MPs
            disp(['Preparing FD measures for: ' stre_txt]);
            temp_txt_fn = fullfile(d, [f '.txt']);
            FD_measures = fmrwhy_qc_calculateFD(MP, options.r, options.FD_threshold); % FD_measures.FD_outliers_regr
            FD = FD_measures.FD;
            if options.FD_threshold == 0
                FD_outliers_regr02 = FD_measures.FD_outliers_regr02;
                FD_outliers_regr05 = FD_measures.FD_outliers_regr05;
                FD_outliers_regr = [FD_outliers_regr02 FD_outliers_regr05];
                col_names_fd = {'framewise_displacement', 'framewise_displacement_censor02', 'framewise_displacement_censor05'};
            else
                FD_outliers_regr = FD_measures.FD_outliers_regr;
                col_names_fd = {'framewise_displacement', 'framewise_displacement_censor'};
            end
            data_table = array2table([FD FD_outliers_regr], 'VariableNames', col_names_fd);
            writetable(data_table, temp_txt_fn, 'Delimiter', '\t');
            [status, msg, msgID] = movefile(temp_txt_fn, options.framewise_displacement_fn);
            disp('Complete!');
            disp('---');
        else
            disp(['FD measures already calculated: ' options.framewise_displacement_fn]);
            disp('---');
            FD_struct = tdfread(options.framewise_displacement_fn);
            FD_measures = struct2array(FD_struct);
            FD = FD_measures(:, 1);
            FD_outliers_regr = FD_measures(:, 2:end);
            col_names_fd = (fieldnames(FD_struct))';
        end
    end

    % -------
    % STEP 2: Nuisance signals from tissue compartment and whole brain masks
    % TODO: figure out from which preprocessed timeseries this has to be calculated. For now using rafunctional.
    % -------
    if include_tissue
        [d, f, e] = fileparts(options.tissue_regr_fn);
        if ~exist(options.tissue_regr_fn, 'file')
            disp(['Calculating tissue compartment signal averages: ' options.tissue_regr_fn]);
            disp('---');
            % Get anatomical tissue compartment masks in individual functional space
            masks = fmrwhy_util_loadMasks(bids_dir, sub);
            % Get timeseries data
            nii = nii_tool('load', options.rafunctional_fn);
            rafunctional_4D = nii.img;
            [Ni, Nj, Nk, Nt] = size(rafunctional_4D);
            rafunctional_2D = reshape(rafunctional_4D, Nj * Nj * Nk, Nt); % [Nj*Nj*Nk, Nt]
            % Get spatial mean per timepoint per tissue compartment
            GM_timeseries = (mean(rafunctional_2D(masks.GM_mask_I, :)))';
            WM_timeseries = (mean(rafunctional_2D(masks.WM_mask_I, :)))';
            CSF_timeseries = (mean(rafunctional_2D(masks.CSF_mask_I, :)))';
            brain_timeseries = (mean(rafunctional_2D(masks.brain_mask_I, :)))';
            % Save to tsv file
            temp_txt_fn = fullfile(d, [f '.txt']);
            col_names = {'grey_matter', 'white_matter', 'csf', 'global_signal'};
            data = [GM_timeseries WM_timeseries CSF_timeseries brain_timeseries];
            data_table = array2table(data, 'VariableNames', col_names);
            writetable(data_table, temp_txt_fn, 'Delimiter', '\t');
            [status, msg, msgID] = movefile(temp_txt_fn, options.tissue_regr_fn);
        else
            disp(['Tissue compartment signal averages already calculated: ' options.tissue_regr_fn]);
            disp('---');
            tissue_regressors = struct2array(tdfread(options.tissue_regr_fn));
            GM_timeseries = tissue_regressors(:, 1);
            WM_timeseries = tissue_regressors(:, 2);
            CSF_timeseries = tissue_regressors(:, 3);
            brain_timeseries = tissue_regressors(:, 4);
        end
    end

    % -------
    % STEP 3: Physiological regressors (Retroicor, HRV, RVT)
    % -------
    if include_physio
        % Grab physio log-file, e.g. sub-001_task-rest_run-1_physio.tsv.gz
        log_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_physio.tsv.gz']);

        % Create struct with options for PhysIO
        physio_options = options.physio.options;
        physio_options.save_dir = fullfile(options.sub_dir_qc, 'func', ['PhysIO_task-' task '_run-' run]);
        physio_options.cardiac_fn = log_fn;
        physio_options.respiration_fn = log_fn;
        physio_options.level = 0; % verbose.level = 0 ==> do not generate figure outputs during batch process

        % Run PhysIO if required
        [d, f, e] = fileparts(options.physio_regr_fn);
        if ~exist(options.physio_regr_fn, 'file')
            disp('Physio regressors not calculated yet. Calculating now.');
            % Run batch process without generating figures (batch calls `tapas_physio_main_create_regressors`)
            phys_data = fmrwhy_batch_PhysIO(physio_options);
            temp_txt_fn = fullfile(d, [f '.txt']);
            col_names = {'retroicor_c1', 'retroicor_c2', 'retroicor_c3', 'retroicor_c4', 'retroicor_c5', 'retroicor_c6', 'retroicor_r1', 'retroicor_r2', 'retroicor_r3', 'retroicor_r4', 'retroicor_r5', 'retroicor_r6', 'retroicor_r7', 'retroicor_r8', 'retroicor_cxr1', 'retroicor_cxr2', 'retroicor_cxr3', 'retroicor_cxr4', 'hrv', 'rvt'};
            data = load(phys_data.multiple_regressors_fn);
            data_table = array2table(data, 'VariableNames', col_names);
            writetable(data_table, temp_txt_fn, 'Delimiter', '\t');
            [status, msg, msgID] = movefile(temp_txt_fn, options.physio_regr_fn);
            % Generate figures (by calling `tapas_physio_review` with updated physio options)
            phys_file = fullfile(physio_options.save_dir, 'physio.mat');
            physmat = load(phys_file);
            physmat.physio.verbose.level = 2;
            physmat.physio.verbose.fig_output_file = ['sub-' sub '_task-' task '_run-' run '_physioQC.jpg'];
            physmat.physio.verbose.show_figs = false;
            physmat.physio.verbose.save_figs = true;
            physmat.physio.verbose.close_figs = true;
            % Run tapas_physio_review
            phys_data_figures = tapas_physio_review(physmat.physio);

        else
            disp('Physio regressors previously calculated. Loading now.');
            phys_data.physio_regr = struct2array(tdfread(options.physio_regr_fn));
        end
    end

    % -------
    % STEP 4: Add all regressors to one array and save to tsv
    % TODO: this does not replace the txt file, i.e. duplicate text file
    % remains. Fix this.
    % -------
    [d, f, e] = fileparts(options.confounds_fn);
    if ~exist(options.confounds_fn, 'file')
        disp('Main confound regressor file not created yet. Creating now.');
        confounds_txt = fullfile(d, [f '.txt']);
        col_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
        data = MP;
        if include_fd
            col_names = [col_names, col_names_fd];
            data = [data FD FD_outliers_regr];
        end
        if include_tissue
            col_names = [col_names, {'white_matter', 'csf', 'global_signal'}];
            data = [data WM_timeseries CSF_timeseries brain_timeseries];
        end
        if include_volterra
            col_names = [col_names, col_names_volterra];
            data = [data MP_volterra];
        end
        if include_physio
            col_names = [col_names, {'retroicor_c1', 'retroicor_c2', 'retroicor_c3', 'retroicor_c4', 'retroicor_c5', 'retroicor_c6', 'retroicor_r1', 'retroicor_r2', 'retroicor_r3', 'retroicor_r4', 'retroicor_r5', 'retroicor_r6', 'retroicor_r7', 'retroicor_r8', 'retroicor_cxr1', 'retroicor_cxr2', 'retroicor_cxr3', 'retroicor_cxr4', 'hrv', 'rvt'}];
            data = [data phys_data.physio_regr];
        end
        dlmwrite(confounds_txt, data, 'delimiter', '\t', 'precision', '%1.7e');
        confounds_tsv = fmrwhy_util_saveAsTSV(confounds_txt, col_names);
    else
        disp('Main confound regressor file already exists.');
    end
