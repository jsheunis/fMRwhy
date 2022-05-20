function confounds_tsv = fmrwhy_bids_preprocMultRegr(bids_dir, sub, task, options, varargin)
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

    % -------------
    % Parse inputs
    % -------------
    filetypes = {'func'};
    descriptions = {'Subject', 'Session', 'Task', 'Acquisition', 'Contrast Enhancing Agent', 'Reconstruction', 'Phase-Encoding Direction', 'Run', 'Echo'};
    entities = {'ses', 'acq', 'rec', 'run', 'echo'}; % these entities are required/optional for func bold data specifically (not other types!)
    formats = {'label', 'label', 'label', 'label', 'label', 'label', 'label', 'index', 'index'};

    validChar = @(x) ischar(x);
    validType = @(x) any(validatestring(x, filetypes));

    p = inputParser;
    addRequired(p, 'bids_dir', validChar);
    addRequired(p, 'sub', validChar);
    addRequired(p, 'task', validChar);
    addRequired(p, 'options');
    for i = 1:numel(entities)
        addParameter(p, entities{i}, '', validChar);
    end
    parse(p, bids_dir, sub, task, options, varargin{:});
    params = p.Results;
    bids_dir = params.bids_dir;
    sub = params.sub;
    task = params.task;
    options = params.options;

    % Default steps to include/exclude
    include_volterra = options.confounds.include_volterra;
    include_fd = options.confounds.include_fd;
    include_tissue = options.confounds.include_tissue;
    include_physio = options.confounds.include_physio;

    % Acces template functional data
    % assignin('base', 'options', options)
    % assignin('base', 'options', params)
    functional_data = bids.query(options.bids_dataset, 'data', 'type', 'bold', 'sub', sub, 'task', task, 'ses', params.ses, 'acq', params.acq, 'rec', params.rec, 'run', params.run);
    % First check if query returns bold data, if not throw error
    if isempty(functional_data)
        msg = 'No functional bold data found for the specified entities';
        fmrwhy_util_createErrorMsg('fmrwhy_bids_preprocFunc', 'dataNotFound', msg);
    else
        % if there is data, check to see how many echos
        is_multiecho = false;
        N_echoes = numel(functional_data);
        if N_echoes > 1
            is_multiecho = true;
        end
    end
    % Assign correct echo if multiecho
    echo = params.echo;
    if is_multiecho
        echo = options.template_echo;
    end

    % Update functional derivate filenames
    options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'echo', echo, 'acq', params.acq, 'rec', params.rec);

    % -------
    % STEP 1: Generate realignment, framewise displacement and censoring regressors
    % -------
    % First realignment parameters
    % Check if this has already been done by seeing if the tsv file with head movement parameters exist
    [d, f, e] = fileparts(options.motion_fn);
    if ~exist(options.motion_fn, 'file')
        % If it does not exist, estimate MPs
        disp(['Estimating 3D realignment parameters for: ' options.current_functional_filename]);
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
            disp(['Preparing FD measures for: ' options.current_functional_filename]);
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
            masks = fmrwhy_util_loadMasks(bids_dir, sub, options, 'ses', params.ses, 'task', task);
            % Get timeseries data from realigned (and slice time corrected, if done)
            if options.include_stc
                nii = nii_tool('load', options.rafunctional_fn);
            else
                nii = nii_tool('load', options.rfunctional_fn);
            end

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
        if 1 == isempty(params.ses)
            log_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_ses-' params.ses '_task-' task '_run-' params.run '_physio.log']);
        else
            log_fn = fullfile(options.sub_dir_preproc, ['ses-' params.ses], 'func', ['sub-' sub '_ses-' params.ses '_task-' task '_run-' params.run '_physio.log']);
        end

        % Create struct with options for PhysIO
        physio_options = options.physio.options;
        if 1 == isempty(params.ses)
            physio_options.save_dir = fullfile(options.sub_dir_qc, 'func', ['PhysIO_task-' task '_run-' params.run]);
        else
            physio_options.save_dir = fullfile(options.sub_dir_qc, ['ses-' params.ses], 'func', ['PhysIO_task-' task '_run-' params.run]);
        end
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
            if 1 == isempty(params.ses)
                physmat.physio.verbose.fig_output_file = ['sub-' sub '_task-' task '_run-' params.run '_physioQC.jpg'];
            else
                physmat.physio.verbose.fig_output_file = ['sub-' sub 'ses' params.ses '_task-' task '_run-' params.run '_physioQC.jpg'];
            end
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
