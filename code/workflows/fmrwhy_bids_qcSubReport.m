function fmrwhy_bids_qcSubReport(bids_dir, sub)


% ------------------------------------
% STEP 1 -- Load default filenames etc
% ------------------------------------

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% Get template echo
echo = options.template_echo;



% ---------------------------
% STEP 2 -- Set up HTML files
% ---------------------------

% Set up datetime strings for filename and content
dt = datetime('now');
[Y,MO,D,H,MI,S] = datevec(dt);
dt_str = [num2str(Y) num2str(MO) num2str(D) num2str(H) num2str(MI) num2str(round(S))];
t = datestr(dt);

% Set up filenames
fmrwhy_dir = options.fmrwhy_dir;
html_template_fn = fullfile(fmrwhy_dir, 'assets', 'fmrwhy_bids_qcSubReportTemplate.htm');
css_template_fn = fullfile(fmrwhy_dir, 'assets', 'fmrwhy_bids_qcReport.css');
js_template_fn = fullfile(fmrwhy_dir, 'assets', 'fmrwhy_bids_qcReport.js');
avatar_template_fn = fullfile(fmrwhy_dir, 'img', 'logo_jsheunis_3.jpeg');

report_dir = fullfile(options.sub_dir_qc, ['report_' dt_str]);
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
% STEP 3 -- Create structure with content to write to HTML report
% ---------------------------------------------------------------
% Locations of assets for html file, e.g. css, js and images.
report.param_asset_css = fullfile('assets', 'fmrwhy_reports.css');
copyfile(css_template_fn, fullfile(report_dir, report.param_asset_css));

report.param_asset_js = fullfile('assets', 'fmrwhy_reports.js');
report.param_avatar = fullfile('assets', 'fmrwhy_logo.jpeg');
copyfile(avatar_template_fn, fullfile(report_dir, report.param_avatar));
% Details about study, subject, data, etc
report.param_sub = sub;
report.param_datetime = t;
% Anatomical montage image locations
brain_mask = fullfile(options.anat_dir_qc, ['sub-' sub '_brain_mask_montage.png']);
gm_mask = fullfile(options.anat_dir_qc, ['sub-' sub '_GM_mask_montage.png']);
wm_mask = fullfile(options.anat_dir_qc, ['sub-' sub '_WM_mask_montage.png']);
csf_mask = fullfile(options.anat_dir_qc, ['sub-' sub '_CSF_mask_montage.png']);
report.param_brain_mask = fullfile('img', ['sub-' sub '_brain_mask_montage.png']);
report.param_gm_mask = fullfile('img', ['sub-' sub '_GM_mask_montage.png']);
report.param_wm_mask = fullfile('img', ['sub-' sub '_WM_mask_montage.png']);
report.param_csf_mask = fullfile('img', ['sub-' sub '_CSF_mask_montage.png']);
copyfile(brain_mask, fullfile(report_dir, report.param_brain_mask));
copyfile(gm_mask, fullfile(report_dir, report.param_gm_mask));
copyfile(wm_mask, fullfile(report_dir, report.param_wm_mask));
copyfile(csf_mask, fullfile(report_dir, report.param_csf_mask));

% ROI montage image locations
% TODO: the ROI filenames are hardcoded in html for now, need to change this in future
roi_img1 = fullfile(options.anat_dir_qc, ['sub-' sub  '_space-individual_desc-leftMotor_roi_montage.png']);
roi_img2 = fullfile(options.anat_dir_qc, ['sub-' sub  '_space-individual_desc-rightMotor_roi_montage.png']);
roi_img3 = fullfile(options.anat_dir_qc, ['sub-' sub  '_space-individual_desc-leftAmygdala_roi_montage.png']);
roi_img4 = fullfile(options.anat_dir_qc, ['sub-' sub  '_space-individual_desc-rightAmygdala_roi_montage.png']);
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

% Stats summary data for specific functional run
tasks_runs = {'rest_run-1', 'motor_run-1', 'emotion_run-1', 'rest_run-2', 'motor_run-2', 'emotion_run-2'};
for i = 1:numel(tasks_runs)
    % write summary metrics
    stats_summary_fn = fullfile(options.func_dir_qc, ['sub-' sub '_task-' tasks_runs{i} '_desc-stats_summary.tsv']);
    stats = tdfread(stats_summary_fn);
    framewise_displacement_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' tasks_runs{i} '_desc-confounds_fd.tsv']);
    fd = tdfread(framewise_displacement_fn);
    report.(['param_fd_mean_' num2str(i)]) = sprintf('%0.2f', mean(fd.framewise_displacement));
    report.(['param_fd_total_' num2str(i)]) = sprintf('%0.2f', sum(fd.framewise_displacement));
    report.(['param_fd_outliers02_' num2str(i)]) = sprintf('%0.0f', numel(find(fd.framewise_displacement_censor02)));
    report.(['param_fd_outliers05_' num2str(i)]) = sprintf('%0.0f', numel(find(fd.framewise_displacement_censor05)));
    report.(['param_zscore_mean_' num2str(i)]) = sprintf('%0.2f', stats.zscore_mean);
    report.(['param_gcor_' num2str(i)]) = sprintf('%0.2f', stats.gcor);
    report.(['param_tsnr_GM_mean_' num2str(i)]) = sprintf('%0.2f', stats.tSNR_mean_GM);
    report.(['param_tsnr_WM_mean_' num2str(i)]) = sprintf('%0.2f', stats.tSNR_mean_WM);
    report.(['param_tsnr_CSF_mean_' num2str(i)]) = sprintf('%0.2f', stats.tSNR_mean_CSF);
    report.(['param_tsnr_brain_mean_' num2str(i)]) = sprintf('%0.2f', stats.tSNR_mean_brain);

    % copy functional qc images
    to_copy.mean_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' tasks_runs{i} '_space-individual_mean.png']);
    to_copy.tsnr_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' tasks_runs{i} '_space-individual_tsnr.png']);
    to_copy.var_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' tasks_runs{i} '_space-individual_var.png']);
    to_copy.std_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' tasks_runs{i} '_space-individual_std.png']);
    % Functional timeseries image locations
    to_copy.rograyplot_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' tasks_runs{i} '_echo-' echo '_desc-RO_grayplot.png']);
    to_copy.gsograyplot_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' tasks_runs{i} '_echo-' echo '_desc-GSO_grayplot.png']);
    % TODO: this is hardcoded for this study, need to generalise it
    to_copy.leftMotorgrayplot_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' tasks_runs{i} '_echo-' echo '_desc-leftMotor_grayplot.png']);
    to_copy.rightMotorgrayplot_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' tasks_runs{i} '_echo-' echo '_desc-rightMotor_grayplot.png']);
    to_copy.leftAmygrayplot_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' tasks_runs{i} '_echo-' echo '_desc-leftAmygdala_grayplot.png']);
    to_copy.rightAmygrayplot_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' tasks_runs{i} '_echo-' echo '_desc-rightAmygdala_grayplot.png']);
    to_copy.biAmygrayplot_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' tasks_runs{i} '_echo-' echo '_desc-bilateralAmygdala_grayplot.png']);
    % PhysIO QC image locations
    to_copy.physioqc_img = fullfile(options.func_dir_qc, ['PhysIO_task-' tasks_runs{i}], ['sub-' sub '_task-' tasks_runs{i} '_physioQC_03.jpg']);
    % Copy all images if they exist
    cp_fields = fieldnames(to_copy);
    for x = 1:numel(cp_fields)
       if exist(to_copy.(cp_fields{x}), 'file')
           copyfile(to_copy.(cp_fields{x}), report_img_dir);
       else
           disp(['File does not exist: ' to_copy.(cp_fields{x})])
       end
    end

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
        fprintf(fid_tmp,'%s\n', newString);
%        fprintf(fid_tmp,'%s', newString);
    end
end
fclose(fid_template);
fclose(fid_tmp);

[status, msg, msgID] = movefile(html_tmp_fn, html_new_fn);
delete(html_tmp_fn)


% -----------------------------------------------------
% STEP 5 -- Read, replace and write data to js file
% -----------------------------------------------------

%

js_strings.param_str2 = [filesep 'sub-' sub '_task-'];
js_strings.param_str3 = '_physioQC_03.jpg';

fields = fieldnames(js_strings);

fid_template = fopen(js_template_fn);
fid_tmp = fopen(js_tmp_fn, 'w');

tline = fgetl(fid_template);
while ischar(tline)
%    disp(tline)
    tline = fgetl(fid_template);
    newString = tline;
    for i = 1:numel(fields)
        newString = strrep(newString, fields{i}, js_strings.(fields{i}));
    end
    if ischar(tline)
        fprintf(fid_tmp,'%s\n', newString);
    end
end
fclose(fid_template);
fclose(fid_tmp);

[status, msg, msgID] = movefile(js_tmp_fn, fullfile(report_dir, report.param_asset_js));
delete(js_tmp_fn)























