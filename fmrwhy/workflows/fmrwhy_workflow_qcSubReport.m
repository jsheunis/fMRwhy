function [report, js_string] = fmrwhy_workflow_qcSubReport(sub, options)
    % A workflow that generates a quality report (in HTML format) for a single subject in a BIDS dataset.
    % This workflow assumes that all necessary preprocessing steps have been completed,
    % including anatomical to functional coregistration, segmentation, ROI registration (if specified),
    % functional time series preprocessing, and the generation of quality control metrics and visualizations.
    % All of these preprocessing steps are run as part of the :func:`fmrwhy_workflow_qc` workflow.
    %
    % :param sub: the subject identifier for which the report is to be generated
    % :param options: the updated `options` structure
    % :type sub: character array
    % :type options: struct
    % :returns: ``[report, js_string]`` - variable array with two strings
    %
    % .. seealso:: This workflow is called by :func:`fmrwhy_workflow_qc`

    % ------------------------------------
    % STEP 0 -- Get BIDS info
    % ------------------------------------
    qc_report_runs = {};
    qc_report_sestasks = {};
    BIDS = bids.layout(options.bids_dir);
    sessions = bids.query(BIDS, 'sessions', 'sub', sub);
    if ~isempty(sessions)
        for s = 1:numel(sessions)
            ses = sessions{s};
            tasks = bids.query(BIDS, 'tasks', 'sub', sub, 'ses', ses);
            % tasks should not be empty
            for t = 1:numel(tasks)
                task = tasks{t};
                qc_report_sestasks = [qc_report_sestasks {['ses-' ses '_task-' task]}];
                runs = bids.query(BIDS, 'runs', 'sub', sub, 'ses', ses, 'task', task);
                for r = 1:numel(runs)
                    rn = runs{r};
                    run_name = ['ses-' ses '_task-' task '_run-' rn];
                    qc_report_runs = [qc_report_runs {run_name}];
                end
            end
        end
    else
        tasks = bids.query(BIDS, 'tasks', 'sub', sub);
        % tasks should not be empty
        for t = 1:numel(tasks)
            task = tasks{t};
            runs = bids.query(BIDS, 'runs', 'sub', sub, 'ses', ses, 'task', task);
            for r = 1:numel(runs)
                rn = runs{r};
                run_name = ['task-' task '_run-' rn];
                qc_report_runs = [qc_report_runs {run_name}];
            end
        end
    end
    options.qc_report_runs = qc_report_runs;
    options.qc_report_sestasks = qc_report_sestasks;

    % ------------------------------------
    % STEP 1 -- Load default filenames etc
    % ------------------------------------

    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_bids_getAnatDerivs(options.bids_dir, sub, options);

    % ---------------------------
    % STEP 2 -- Set up HTML, CSS, JS files
    % ---------------------------

    % Set up datetime strings for filename and content
    dt = datetime('now');
    [Y, MO, D, H, MI, S] = datevec(dt);
    dt_str = [num2str(Y) sprintf('%02d', MO) sprintf('%02d', D) sprintf('%02d', H) sprintf('%02d', MI) sprintf('%02d', round(S))];
    t = datestr(dt);

    % Set up filenames
    fmrwhy_dir = options.fmrwhy_dir;
    html_template_fn = fullfile(fmrwhy_dir, 'assets', 'fmrwhy_bids_qcSubReportTemplate.htm');
    css_template_fn = fullfile(fmrwhy_dir, 'assets', 'fmrwhy_bids_qcReport.css');
    js_template_fn = fullfile(fmrwhy_dir, 'assets', 'fmrwhy_bids_qcSubReportTemplate.js');
    avatar_template_fn = fullfile(fmrwhy_dir, 'img', 'logo_jsheunis_3.jpeg');

    report_dir = fullfile(options.qc_dir, ['sub-' sub], ['report_' dt_str]);
    report_assets_dir = fullfile(report_dir, 'assets');
    report_img_dir = fullfile(report_dir, 'img');
    if ~exist(report_dir, 'dir')
        mkdir(report_dir);
    end
    if ~exist(report_assets_dir, 'dir')
        mkdir(report_assets_dir);
    end
    if ~exist(report_img_dir, 'dir')
        mkdir(report_img_dir);
    end

    html_new_fn = fullfile(report_dir, ['sub-' sub '_desc-QCreport_' dt_str '.html']);
    html_tmp_fn = 'temp.html';
    js_tmp_fn = 'temp.js';

    % ---------------------------------------------------------------
    % STEP 3 -- Create structures with content to write to HTML report
    % ---------------------------------------------------------------
    report = struct;
    bids_dataset = struct;

    % Locations of assets for html file, e.g. css, js and images.
    report.param_asset_js = fullfile('assets', 'fmrwhy_bids_qcReport.js'); % don't copy js file yet; first to be updated.
    report.param_asset_css = fullfile('assets', 'fmrwhy_bids_qcReport.css');
    copyfile(css_template_fn, fullfile(report_dir, report.param_asset_css));
    report.param_avatar = fullfile('assets', 'fmrwhy_logo.jpeg');
    copyfile(avatar_template_fn, fullfile(report_dir, report.param_avatar));

    % Details about study, subject, data, etc
    report.param_dataset_name = options.qc_dataset_name;
    report.param_sub = sub;
    report.param_datetime = t;
    report.param_anat_res = options.qc_anat_res;
    report.param_func_res = options.qc_func_res;
    report.param_func_acq = options.qc_func_acq;
    report.param_func_runs = options.qc_func_runs;
    report.param_first_run = options.qc_report_runs{1};
    report.param_first_sestask = '';
    if ~isempty(options.qc_report_sestasks)
        report.param_first_sestask = options.qc_report_sestasks{1};
    end

    % Populate bids_dataset structure with existing variables
    bids_dataset.sub = ['sub-' sub];
    bids_dataset.sessions = options.sessions;
    bids_dataset.tasks = options.tasks;
    bids_dataset.runs = options.runs;
    bids_dataset.anat_template_session = options.anat_template_session;
    bids_dataset.template_session = options.template_session;
    bids_dataset.template_task = options.template_task;
    bids_dataset.template_run = options.template_run;
    bids_dataset.template_echo = options.template_echo;
    bids_dataset.has_sessions = options.has_sessions;
    bids_dataset.has_runs = options.has_runs;
    bids_dataset.is_multiecho = options.is_multiecho;
    bids_dataset.map_rois = false;
    if options.map_rois
        bids_dataset.map_rois = true;
    end
    bids_dataset.include_physio = false;
    if options.confounds.include_physio
        bids_dataset.include_physio = true;
    end
    bids_dataset.physio_str = '_physioQC_03.jpg';

    % Copy all PNG files recursively to report image directory
    sub_qc_dir = fullfile(options.qc_dir, ['sub-' sub]);
    png_list = dir(fullfile(sub_qc_dir, '**/*.png'));  % get list of png files and folders in all subfolder
    png_list = png_list(~[png_list.isdir]);  % remove folders from list

    for i = 1:numel(png_list)
        % If sessions exist
        fn = fullfile(png_list(i).folder, png_list(i).name);
        to_fn = fullfile(report_img_dir, png_list(i).name);
        copyfile(fn, to_fn);
    end

    % [filename, filepath] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'ext', '_T1w.nii');
    % Anatomical montage image locations - all anatomical QC outputs should be located in the 'anat' dir (in QC derivatives) of the template session; if no sessions ==> template session is the main 'anat' dir
    % brain_mask = fullfile(options.qc_dir, filepath, ['sub-' sub '_brain_mask_montage.png']);
    % gm_mask = fullfile(options.qc_dir, filepath, ['sub-' sub '_GM_mask_montage.png']);
    % wm_mask = fullfile(options.qc_dir, filepath, ['sub-' sub '_WM_mask_montage.png']);
    % csf_mask = fullfile(options.qc_dir, filepath, ['sub-' sub '_CSF_mask_montage.png']);
    % report.param_brain_mask = fullfile('img', ['sub-' sub '_space-individual_desc-brain_mask_montage.png']);
    % report.param_gm_mask = fullfile('img', ['sub-' sub '_space-individual_desc-GM_mask_montage.png']);
    % report.param_wm_mask = fullfile('img', ['sub-' sub '_space-individual_desc-WM_mask_montage.png']);
    % report.param_csf_mask = fullfile('img', ['sub-' sub '_space-individual_desc-CSF_mask_montage.png']);

    % TODO: this is the intended code for the correct version of the report, the code above is used for rt-me-fMRI dataset (beta version, v0.0.1)
    % brain_mask = fullfile(options.qc_dir, filepath, ['sub-' sub '_space-individual_desc-brain_mask_montage.png']);
    % gm_mask = fullfile(options.qc_dir, filepath, ['sub-' sub '_space-individual_desc-GM_mask_montage.png']);
    % wm_mask = fullfile(options.qc_dir, filepath, ['sub-' sub '_space-individual_desc-WM_mask_montage.png']);
    % csf_mask = fullfile(options.qc_dir, filepath, ['sub-' sub '_space-individual_desc-CSF_mask_montage.png']);
    % report.param_brain_mask = fullfile('img', ['sub-' sub '_space-individual_desc-brain_mask_montage.png']);
    % report.param_gm_mask = fullfile('img', ['sub-' sub '_space-individual_desc-GM_mask_montage.png']);
    % report.param_wm_mask = fullfile('img', ['sub-' sub '_space-individual_desc-WM_mask_montage.png']);
    % report.param_csf_mask = fullfile('img', ['sub-' sub '_space-individual_desc-CSF_mask_montage.png']);
    % copyfile(brain_mask, fullfile(report_dir, report.param_brain_mask));
    % copyfile(gm_mask, fullfile(report_dir, report.param_gm_mask));
    % copyfile(wm_mask, fullfile(report_dir, report.param_wm_mask));
    % copyfile(csf_mask, fullfile(report_dir, report.param_csf_mask));

    % ROI montage image locations - not working yet
    if options.map_rois
        roi_img1 = fullfile(options.qc_dir, filepath, ['sub-' sub  '_space-individual_desc-leftMotor_roi_montage.png']);
        roi_img2 = fullfile(options.qc_dir, filepath, ['sub-' sub  '_space-individual_desc-rightMotor_roi_montage.png']);
        roi_img3 = fullfile(options.qc_dir, filepath, ['sub-' sub  '_space-individual_desc-leftAmygdala_roi_montage.png']);
        roi_img4 = fullfile(options.qc_dir, filepath, ['sub-' sub  '_space-individual_desc-rightAmygdala_roi_montage.png']);
        report.param_roi_img1 = fullfile('img', ['sub-' sub  '_space-individual_desc-leftMotor_roi_montage.png']);
        report.param_roi_img2 = fullfile('img', ['sub-' sub  '_space-individual_desc-rightMotor_roi_montage.png']);
        report.param_roi_img3 = fullfile('img', ['sub-' sub  '_space-individual_desc-leftAmygdala_roi_montage.png']);
        report.param_roi_img4 = fullfile('img', ['sub-' sub  '_space-individual_desc-rightAmygdala_roi_montage.png']);
        copyfile(roi_img1, fullfile(report_dir, report.param_roi_img1));
        copyfile(roi_img2, fullfile(report_dir, report.param_roi_img2));
        copyfile(roi_img3, fullfile(report_dir, report.param_roi_img3));
        copyfile(roi_img4, fullfile(report_dir, report.param_roi_img4));
        report.param_roi_name1 = 'Left Motor';
        report.param_roi_name2 = 'Right Motor';
        report.param_roi_name3 = 'Left Amygdala';
        report.param_roi_name4 = 'Right Amygdala';

        % For report purposes, which ROIs to include in grayplots. TODO.
        bids_dataset.roi_desc = {'leftMotor', 'biAmygdala'};
    end

    % Stats summary data for each specific functional run, specified in settings file
    bids_dataset.report_runs = options.qc_report_runs;
    bids_dataset.report_sestasks = options.qc_report_sestasks;
    tasks_runs = options.qc_report_runs;
    bids_dataset.all_runs_stats = cell(1, numel(tasks_runs));
    % [func_filename, func_filepath] = fmrwhy_bids_constructFilename('func', 'sub', sub, 'ses', 'task', options.template_task, 'run', options.template_run, 'space', 'individual', 'ext', '_bold.nii');
    for i = 1:numel(tasks_runs)
        if contains(tasks_runs{i}, 'ses-')
            ses_arr = extractBetween(tasks_runs{i}, 'ses-', '_');
            ses = ses_arr{1};
            f_path_qc = fullfile(options.qc_dir, ['sub-' sub], ['ses-' ses], 'func');
            f_path_preproc = fullfile(options.preproc_dir, ['sub-' sub], ['ses-' ses], 'func');
        else
            f_path_qc = fullfile(options.qc_dir, ['sub-' sub], 'func');
            f_path_preproc = fullfile(options.preproc_dir, ['sub-' sub], 'func');

        end

        % write summary metrics
        stats_summary_fn = fullfile(f_path_qc, ['sub-' sub '_' tasks_runs{i} '_desc-stats_summary.tsv']);
        stats = tdfread(stats_summary_fn);
        framewise_displacement_fn = fullfile(f_path_preproc, ['sub-' sub '_' tasks_runs{i} '_desc-confounds_fd.tsv']);
        fd = tdfread(framewise_displacement_fn);

        bids_dataset.all_runs_stats{i}.mean_fd = sprintf('%0.2f', mean(fd.framewise_displacement));
        bids_dataset.all_runs_stats{i}.total_fd = sprintf('%0.2f', sum(fd.framewise_displacement));
        bids_dataset.all_runs_stats{i}.fd_outliers02 = sprintf('%0.0f', numel(find(fd.framewise_displacement_censor02)));
        bids_dataset.all_runs_stats{i}.fd_outliers05 = sprintf('%0.0f', numel(find(fd.framewise_displacement_censor05)));
        bids_dataset.all_runs_stats{i}.mean_zscore = sprintf('%0.2f', stats.zscore_mean);
        bids_dataset.all_runs_stats{i}.gcor = sprintf('%0.2f', stats.gcor);
        bids_dataset.all_runs_stats{i}.tsnr_gm = sprintf('%0.2f', stats.tSNR_mean_GM);
        bids_dataset.all_runs_stats{i}.tsnr_wm = sprintf('%0.2f', stats.tSNR_mean_WM);
        bids_dataset.all_runs_stats{i}.tsnr_csf = sprintf('%0.2f', stats.tSNR_mean_CSF);
        bids_dataset.all_runs_stats{i}.tsnr_brain = sprintf('%0.2f', stats.tSNR_mean_brain);

        % % copy functional qc images
        % to_copy.mean_img = fullfile(options.qc_dir, func_filepath, ['sub-' sub '_' tasks_runs{i} '_space-individual_mean.png']);
        % to_copy.tsnr_img = fullfile(options.qc_dir, func_filepath, ['sub-' sub '_' tasks_runs{i} '_space-individual_tsnr.png']);
        % to_copy.var_img = fullfile(options.qc_dir, func_filepath, ['sub-' sub '_' tasks_runs{i} '_space-individual_var.png']);
        % to_copy.std_img = fullfile(options.qc_dir, func_filepath, ['sub-' sub '_' tasks_runs{i} '_space-individual_std.png']);
        % % Functional timeseries image locations
        % if options.is_multiecho
        %     to_copy.rograyplot_img = fullfile(options.qc_dir, func_filepath, ['sub-' sub '_' tasks_runs{i} '_echo-' options.template_echo '_desc-RO_grayplot.png']);
        %     to_copy.gsograyplot_img = fullfile(options.qc_dir, func_filepath, ['sub-' sub '_' tasks_runs{i} '_echo-' options.template_echo '_desc-GSO_grayplot.png']);
        % else
        %     to_copy.rograyplot_img = fullfile(options.qc_dir, func_filepath, ['sub-' sub '_' tasks_runs{i} '_desc-RO_grayplot.png']);
        %     to_copy.gsograyplot_img = fullfile(options.qc_dir, func_filepath, ['sub-' sub '_' tasks_runs{i} '_desc-GSO_grayplot.png']);
        % end

        % % TODO: this is hardcoded for this study, need to generalise it
        % if options.map_rois
        %     to_copy.leftMotorgrayplot_img = fullfile(options.qc_dir, func_filepath, ['sub-' sub '_' tasks_runs{i} '_echo-' options.template_echo '_desc-leftMotor_grayplot.png']);
        %     to_copy.rightMotorgrayplot_img = fullfile(options.qc_dir, func_filepath, ['sub-' sub '_' tasks_runs{i} '_echo-' options.template_echo '_desc-rightMotor_grayplot.png']);
        %     to_copy.leftAmygrayplot_img = fullfile(options.qc_dir, func_filepath, ['sub-' sub '_' tasks_runs{i} '_echo-' options.template_echo '_desc-leftAmygdala_grayplot.png']);
        %     to_copy.rightAmygrayplot_img = fullfile(options.qc_dir, func_filepath, ['sub-' sub '_' tasks_runs{i} '_echo-' options.template_echo '_desc-rightAmygdala_grayplot.png']);
        %     to_copy.biAmygrayplot_img = fullfile(options.qc_dir, func_filepath, ['sub-' sub '_' tasks_runs{i} '_echo-' options.template_echo '_desc-bilateralAmygdala_grayplot.png']);
        % end

        % if options.confounds.include_physio
        %     % PhysIO QC image locations
        %     to_copy.physioqc_img = fullfile(options.qc_dir, func_filepath, ['PhysIO_' tasks_runs{i}], ['sub-' sub '_' tasks_runs{i} '_physioQC_03.jpg']);
        % end

        % Copy all images if they exist
        % cp_fields = fieldnames(to_copy);
        % for x = 1:numel(cp_fields)
        %     if exist(to_copy.(cp_fields{x}), 'file')
        %         copyfile(to_copy.(cp_fields{x}), report_img_dir);
        %     else
        %         disp(['File does not exist: ' to_copy.(cp_fields{x})]);
        %     end
        % end

    end

    % -----------------------------------------------------
    % STEP 4 -- Read, replace and write data to HTML report
    % -----------------------------------------------------

    %
    fields = fieldnames(report);

    fid_template = fopen(html_template_fn);
    fid_tmp = fopen(html_tmp_fn, 'w');

    tline = fgetl(fid_template);
    while ischar(tline)
        %    disp(tline)
        tline = fgetl(fid_template);
        newString = tline;
        for i = 1:numel(fields)
            newString = strrep(newString, fields{i}, report.(fields{i}));
        end
        if ischar(tline)
            fprintf(fid_tmp, '%s\n', newString);
            %        fprintf(fid_tmp,'%s', newString);
        end
    end
    fclose(fid_template);
    fclose(fid_tmp);

    [status, msg, msgID] = movefile(html_tmp_fn, html_new_fn);
    delete(html_tmp_fn);

    % -----------------------------------------------------
    % STEP 5 -- Read, replace and write data to js file
    % -----------------------------------------------------

    js_string = jsonencode(bids_dataset);

    fid_template = fopen(js_template_fn);
    fid_tmp = fopen(js_tmp_fn, 'w');

    tline = fgetl(fid_template);
    while ischar(tline)
        %    disp(tline)
        tline = fgetl(fid_template);
        newString = tline;
        newString = strrep(newString, 'param_bids_dataset', ['bids_dataset = ' js_string]);
        if ischar(tline)
            fprintf(fid_tmp, '%s\n', newString);
        end
    end
    fclose(fid_template);
    fclose(fid_tmp);

    [status, msg, msgID] = movefile(js_tmp_fn, fullfile(report_dir, report.param_asset_js));
    delete(js_tmp_fn);
