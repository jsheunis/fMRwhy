function fmrwhy_bids_preprocFunc(bids_dir, sub, task, options, varargin)
    % --------------------------------------------------------------------------

    % Copyright statement....

    % --------------------------------------------------------------------------
    % DEFINITION
    % --------------------------------------------------------------------------
    % Function to run basic functional preprocessing steps that are required for
    % several subsequent analysis steps and quality control.

    % THIS PIPELINE IS RUN ON A SIGNLE FUNCTIONAL TIMESERIES (could be single or multi-echo)

    % STEPS:

    % QUESTION: should functional localisers also be done based on combined echo data? Perhaps this is worth another research question?

    % INPUT:

    % OUTPUT:

    % --------------------------------------------------------------------------

    % -------------
    % Parse inputs
    % -------------
    filetypes = {'func'};
    descriptions = {'Session', 'Acquisition', 'Contrast Enhancing Agent', 'Reconstruction', 'Phase-Encoding Direction', 'Run', 'Echo'};
    entities = {'ses', 'acq', 'rec', 'run', 'echo'}; % these entities are required/optional for func bold data specifically (not other types!)
    formats = {'label', 'label', 'label', 'label', 'label', 'index', 'index'};

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

    % -------------
    % Run fmrwhy_bids_preprocFunc
    % -------------
    disp('---');
    disp('*** Running fmrwhy_bids_preprocFunc ***');
    disp('---');
    disp('---');

    %% Update workflow params with subject functional derivative filenames
    % options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, ses, task, run, echo, options);

    % Grab relevant functional data
    functional_data = bids.query(options.bids_dataset, 'data', 'type', 'bold', 'sub', params.sub, 'ses', params.ses, 'task', params.task, 'acq', params.acq, 'rec', params.rec, 'run', params.run);
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

    % -------
    % STEP 1: Estimate 3D volume realignment parameters from raw data
    % -------
    % First access template timeseries information
    echo = params.echo;
    if is_multiecho
        echo = options.template_echo;
    end

    % Update functional derivate filenames
    options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'echo', echo, 'acq', params.acq, 'rec', params.rec);

    % Check if realignment has already been done by seeing if the tsv file with head movement parameters exist
    [d, f, e] = fileparts(options.motion_fn);
    if ~exist(options.motion_fn, 'file')
        % If it does not exist estimate MPs
        disp(['Estimating 3D realignment parameters for: ' options.current_functional_filename]);
        realign_measures = fmrwhy_batch_realignEst(options.functional_fn, options.template_fn);
        temp_txt_fn = fullfile(d, [f '.txt']);
        col_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
        data = load(realign_measures.mp_fn);
        data_table = array2table(data, 'VariableNames', col_names);
        writetable(data_table, temp_txt_fn, 'Delimiter', '\t');
        [status, msg, msgID] = movefile(temp_txt_fn, options.motion_fn);
        disp('Complete!');
        disp('---');
    else
        disp(['3D realignment parameters already estimated: ' options.motion_fn]);
        disp('---');
    end

    %%
    % -------
    % STEP 2: Slice timing correction (default is exlcude)
    % -------
    if options.include_stc
        if is_multiecho
            for e = 1:numel(N_echoes)
                % Update workflow params with subject functional derivative filenames
                options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'echo', e, 'acq', params.acq, 'rec', params.rec);
                if ~exist(options.afunctional_fn, 'file')
                    disp(['Performing slice timing correction on: ' options.current_functional_filename]);
                    fmrwhy_batch_sliceTiming(options.functional_fn, options.afunctional_fn, options);
                    disp('Complete!');
                    disp('---');
                else
                    disp(['Slice timing correction already completed for: ' options.current_functional_filename]);
                    disp('---');
                end
            end
        else
            % Update workflow params with subject functional derivative filenames
            options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'acq', params.acq, 'rec', params.rec);
            if ~exist(options.afunctional_fn, 'file')
                disp(['Performing slice timing correction on: ' options.current_functional_filename]);
                fmrwhy_batch_sliceTiming(options.functional_fn, options.afunctional_fn, options);
                disp('Complete!');
                disp('---');
            else
                disp(['Slice timing correction already completed for: ' options.current_functional_filename]);
                disp('---');
            end
        end
    end

    %%
    % -------
    % STEP 3: 3D volume realignment (tsnr and other stat measures use slice time corrected and realigned data)
    % -------
    motion_struct = tdfread(options.motion_fn);
    motion_params = struct2array(motion_struct);
    if is_multiecho
        for e = 1:numel(N_echoes)
            % Update workflow params with subject functional derivative filenames
            options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'echo', e, 'acq', params.acq, 'rec', params.rec);
            % Realign raw timeseries data
            if ~exist(options.rfunctional_fn, 'file')
                disp(['Performing 3D realignment on raw timeseries: ' options.current_functional_filename]);
                fmrwhy_util_applyTransform(options.functional_fn, motion_params, options.template_fn, options.rfunctional_fn);
                disp('Complete!');
                disp('---');
            else
                disp(['3D realignment already completed for raw timeseries: ' options.current_functional_filename]);
                disp('---');
            end
            % Realign slice time corrected timeseries data
            if options.include_stc
                if ~exist(options.rafunctional_fn, 'file')
                    disp(['Performing 3D realignment on slice time corrected timeseries: ' options.current_functional_filename]);
                    fmrwhy_util_applyTransform(options.afunctional_fn, motion_params, options.template_fn, options.rafunctional_fn);
                    disp('Complete!');
                    disp('---');
                else
                    disp(['3D realignment already completed for slice time corrected timeseries: ' options.current_functional_filename]);
                    disp('---');
                end
            end
        end
    else
        % Update workflow params with subject functional derivative filenames
        options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'acq', params.acq, 'rec', params.rec);
        % Realign raw timeseries data
        if ~exist(options.rfunctional_fn, 'file')
            disp(['Performing 3D realignment on raw timeseries: ' options.current_functional_filename]);
            fmrwhy_util_applyTransform(options.functional_fn, motion_params, options.template_fn, options.rfunctional_fn);
            disp('Complete!');
            disp('---');
        else
            disp(['3D realignment already completed for raw timeseries: ' options.current_functional_filename]);
            disp('---');
        end
        % Realign slice time corrected timeseries data
        if options.include_stc
            if ~exist(options.rafunctional_fn, 'file')
                disp(['Performing 3D realignment on slice time corrected timeseries: ' options.current_functional_filename]);
                fmrwhy_util_applyTransform(options.afunctional_fn, motion_params, options.template_fn, options.rafunctional_fn);
                disp('Complete!');
                disp('---');
            else
                disp(['3D realignment already completed for slice time corrected timeseries: ' options.current_functional_filename]);
                disp('---');
            end
        end
    end

    %%
    % -------
    % STEP 4: spatial smoothing (the plot uses smoothed data otherwise unprocessed)
    % -------
    if is_multiecho
        for e = 1:numel(N_echoes)
            % Update workflow params with subject functional derivative filenames
            options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'echo', e, 'acq', params.acq, 'rec', params.rec);
            % Smooth raw timeseries data
            if ~exist(options.sfunctional_fn, 'file')
                disp(['Performing spatial smoothing on raw timeseries: ' options.current_functional_filename]);
                fmrwhy_batch_smooth(options.functional_fn, options.sfunctional_fn, options.fwhm);
                disp('Complete!');
                disp('---');
            else
                disp(['Spatial smoothing already completed for raw timeseries: ' options.current_functional_filename]);
                disp('---');
            end
            if basicfunc_full
                % Smooth realigned timeseries data
                if ~exist(options.srfunctional_fn, 'file')
                    disp(['Performing spatial smoothing on realigned timeseries: ' options.current_functional_filename]);
                    fmrwhy_batch_smooth(options.rfunctional_fn, options.srfunctional_fn, options.fwhm);
                    disp('Complete!');
                    disp('---');
                else
                    disp(['Spatial smoothing already completed for realigned timeseries: ' options.current_functional_filename]);
                    disp('---');
                end
                if options.include_stc
                    % Smooth realigned and slice time corrected timeseries data
                    if ~exist(options.srafunctional_fn, 'file')
                        disp(['Performing spatial smoothing on realigned and slice time corrected timeseries: ' options.current_functional_filename]);
                        fmrwhy_batch_smooth(options.rafunctional_fn, options.srafunctional_fn, options.fwhm);
                        disp('Complete!');
                        disp('---');
                    else
                        disp(['Spatial smoothing already completed for realigned and slice time corrected timeseries: ' options.current_functional_filename]);
                        disp('---');
                    end
                end
            end
        end
    else
        % Update workflow params with subject functional derivative filenames
        options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'acq', params.acq, 'rec', params.rec);
        % Smooth raw timeseries data
        if ~exist(options.sfunctional_fn, 'file')
            disp(['Performing spatial smoothing on raw timeseries: ' options.current_functional_filename]);
            fmrwhy_batch_smooth(options.functional_fn, options.sfunctional_fn, options.fwhm);
            disp('Complete!');
            disp('---');
        else
            disp(['Spatial smoothing already completed for raw timeseries: ' options.current_functional_filename]);
            disp('---');
        end
        if options.basicfunc_full
            % Smooth realigned timeseries data
            if ~exist(options.srfunctional_fn, 'file')
                disp(['Performing spatial smoothing on realigned timeseries: ' options.current_functional_filename]);
                fmrwhy_batch_smooth(options.rfunctional_fn, options.srfunctional_fn, options.fwhm);
                disp('Complete!');
                disp('---');
            else
                disp(['Spatial smoothing already completed for realigned timeseries: ' options.current_functional_filename]);
                disp('---');
            end
            if options.include_stc
                % Smooth realigned and slice time corrected timeseries data
                if ~exist(options.srafunctional_fn, 'file')
                    disp(['Performing spatial smoothing on realigned and slice time corrected timeseries: ' options.current_functional_filename]);
                    fmrwhy_batch_smooth(options.rafunctional_fn, options.srafunctional_fn, options.fwhm);
                    disp('Complete!');
                    disp('---');
                else
                    disp(['Spatial smoothing already completed for realigned and slice time corrected timeseries: ' options.current_functional_filename]);
                    disp('---');
                end
            end
        end
    end

    %%
    % -------
    % STEP 5: Generate multiple regressors for GLM analysis and QC. In case of multiecho: all from template echo timeseries
    % Includes: 3D realignment parameters, framewise displacement, FD censoring, tissue compartment signals, retroicor and HRV+RVT
    % -------
    if is_multiecho
        e = options.template_echo;
    else
        e = params.echo;
    end
    options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'echo', e, 'acq', params.acq, 'rec', params.rec);

    if ~exist(options.confounds_fn, 'file')
        disp(['Generating multiple regressors for GLM analysis and QC']);
        fmrwhy_bids_preprocMultRegr(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'echo', e, 'acq', params.acq, 'rec', params.rec);
        disp('Complete!');
        disp('---');
    else
        disp(['Multiple regressors already generated: ' options.confounds_fn]);
        disp('---');
    end
