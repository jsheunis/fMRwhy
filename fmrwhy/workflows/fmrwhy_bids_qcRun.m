function fmrwhy_bids_qcRun(bids_dir, sub, task, options, varargin)

    % Software dependences:
    % - Matlab vX
    % - SPM12 r7771
    % - (TAPAS PhysIO Toolbox vX)

    % Data required in order for this function to run as intended:
    % - T1w (raw data)
    % - Single run of fMRI timeseries (raw data)
    % - Head motion/movement parameters derived from unprocessed data (i.e. from realignment of raw fMRI timeseries)
    % - (T1w coregistered to template functional volume (from fmrwhy_preproc_structFunc.m))
    % - (Segmentations of GM, WM, CSF in template functional volume space (from fmrwhy_preproc_structFunc.m))

    % Code steps:
    % 1. Get BIDS directory of run
    % 2. Create relevant derivative directories (qc and preproc) if it doesn't exist yet
    % 3. QC steps for fmri runs in subject folder:
    %   - Head motion/movement parameters derived from unprocessed data (i.e. from realignment of raw fMRI timeseries)
    %   - Framewise displacement
    %   - Statistical measures / images (tsnr, variance, std, psc, DVARS)
    %   - The plot
    %   - PhysIO quality metrics and plots

    % INPUT
    % - flag to generate pictures or not
    % - flag to generate nii images or not

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

    % ------------------------
    % SECTION A: ANATOMICAL QC
    % ------------------------

    options = fmrwhy_bids_getAnatDerivs(bids_dir, sub, options, 'ses', params.ses, 'task', params.task);

    % -------
    % STEP 1: Contours of tissue masks on mean EPI / template EPI (/ anatomical image?)
    % -------
    % TODO: this should hide figures and only print them to png. currently figures are popping up.
    mask_desc = {'GM', 'WM', 'CSF', 'brain'};
    for i = 1:numel(mask_desc)
        if isempty(options.anat_template_session)
            [filename, filepath] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'ses', params.ses, 'task', params.task, 'space', 'individual', 'desc', mask_desc{i}, 'ext', '_mask_montage.png');
        else
            [filename, filepath] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'ses', options.anat_template_session, 'task', params.task, 'space', 'individual', 'desc', mask_desc{i}, 'ext', '_mask_montage.png');
        end
        mask_montage_fns{i} = fullfile(options.qc_dir, filepath, filename);
    end
    run_montage = 0;
    for i = 1:numel(mask_montage_fns)
        if ~exist(mask_montage_fns{i}, 'file')
            disp(['Mask overlay montage does not exist yet: ' mask_montage_fns{i}]);
            run_montage = 1;
        end
    end
    if run_montage || options.qc_overwrite_tissuecontours
        disp('Creating mask overlay montages');
        montage_data = fmrwhy_bids_qcCreateMaskMontages(bids_dir, sub, 1, options, 'ses', params.ses, 'task', params.task);
        disp('Complete!');
        disp('---');
    else
        disp('All mask overlay montages exist.');
        disp('---');
    end

    % -------
    % STEP 2: Contours of anatomical ROIs on mean EPI / template EPI (/ anatomical image?)
    % -------

    if options.map_rois
        [p1, frm1, rg1, dim1] = fmrwhy_util_readOrientNifti(options.rcoregest_anatomical_fn);
        % Loop through all tasks in BIDS structure
        for i = 1:numel(options.tasks)
            % Ignore the 'rest' task (assume there is no task ROI for this; have to change in future if RSnetworks available to be normalised or something)
            if strcmp(options.tasks{i}, 'rest') ~= 1
                % Loop through all ROIs for the particular task
                for j = 1:numel(options.roi.(options.tasks{i}).orig_fn)
                    [p2, frm2, rg2, dim2] = fmrwhy_util_readOrientNifti(options.roi.(options.tasks{i}).rroi_fn{j});
                    overlay_img = fmrwhy_util_createBinaryImg(p2.nii.img, 0.1);
                    title = options.roi.(options.tasks{i}).name{j};
                    if isempty(options.anat_template_session)
                        [filename, filepath] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'ses', params.ses, 'task', params.task, 'space', 'individual', 'desc', options.roi.(options.tasks{i}).desc{j}, 'ext', '_roi_montage.png');
                    else
                        [filename, filepath] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'ses', options.anat_template_session, 'task', params.task, 'space', 'individual', 'desc', options.roi.(options.tasks{i}).desc{j}, 'ext', '_roi_montage.png');
                    end
                    saveAs_fn = fullfile(options.qc_dir, filepath, filename);
                    if ~exist(saveAs_fn, 'file') || options.qc_overwrite_ROIcontours
                        fmrwhy_util_createOverlayMontage(p1.nii.img, overlay_img, 9, 1, title, 'gray', 'off', 'max', [], [255, 0, 0], saveAs_fn);
                    else
                        disp(['File already exists: ' saveAs_fn]);
                    end
                end
            end
        end
    end

    % ------------------------
    % SECTION B: FUNCTIONAL QC
    % ------------------------

    % -------
    % STEP 1: Grab functionally relevant data, test for multi-echo
    % -------
    % TODO: decide if this needs to be done here, or if the single correct functional run should be sent to this function instead. Going with the latter for now.

    % -------
    % STEP 2: Calculate and generate statistical measures and images (tsnr, variance, std, psc, DVARS)
    % -------
    options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'echo', params.echo, 'acq', params.acq, 'rec', params.rec);

    stats = fmrwhy_bids_qcCreateStatsOutput(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'echo', params.echo, 'acq', params.acq, 'rec', params.rec);

    % -------
    % STEP 3: Create The Plot for whole brain
    % -------
    theplot_fn = {};
    [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', sub, 'ses', params.ses, 'task', task, 'run', params.run, 'echo', params.echo, 'acq', params.acq, 'rec', params.rec, 'desc', 'RO', 'ext', '_grayplot.png');
    theplot_fn{1} = fullfile(options.qc_dir, filepath, filename);
    [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', sub, 'ses', params.ses, 'task', task, 'run', params.run, 'echo', params.echo, 'acq', params.acq, 'rec', params.rec, 'desc', 'GSO', 'ext', '_grayplot.png');
    theplot_fn{2} = fullfile(options.qc_dir, filepath, filename);
    if ~exist(theplot_fn{1}, 'file') || ~exist(theplot_fn{2}, 'file') || options.qc_overwrite_theplot
        options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'echo', params.echo, 'acq', params.acq, 'rec', params.rec);
        fmrwhy_bids_qcCreateThePlot(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'echo', params.echo, 'acq', params.acq, 'rec', params.rec);
    end

    % -------
    % STEP 4: Create The Plot for rois
    % -------
    if options.map_rois
        % Update workflow params with subject functional derivative filenames
        options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, task, options, 'ses', params.ses, 'run', params.run, 'echo', params.echo, 'acq', params.acq, 'rec', params.rec);
        % Ignore the 'rest' task (assume there is no task ROI for this; have to change in future if RSnetworks available to be normalised or something)
        if strcmp(task, 'rest') ~= 1
            % Loop through all ROIs for the particular task
            for j = 1:numel(options.roi.(task).orig_fn)
                functional_fn = options.sfunctional_fn;
                desc = options.roi.(task).desc{j};
                [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', sub, 'ses', params.ses, 'task', task, 'run', params.run, 'echo', params.echo, 'acq', params.acq, 'rec', params.rec, 'desc', desc, 'ext', '_grayplot.png');
                saveAs_fn = fullfile(options.qc_dir, filepath, filename);

                task_info.TR = options.firstlevel.(task).(['run' run]).sess_params.timing_RT;
                task_info.onsets = options.firstlevel.(task).(['run' run]).plot_params.cond_onset;
                task_info.durations = options.firstlevel.(task).(['run' run]).plot_params.cond_duration;
                task_info.precision = 1;

                if ~exist(saveAs_fn, 'file') || options.qc_overwrite_theplot
                    trace_info = [];
                    fmrwhy_util_thePlotROI(functional_fn, options.brain_mask_fn, options.roi.(task).rroi_fn{j}, task_info, trace_info, saveAs_fn);
                else
                    disp(['File already exists: ' saveAs_fn]);
                end
            end
        else
            disp('---');
            disp('Not creating ROI timeseries plots for task = rest.');
            disp('---');
        end
    end
