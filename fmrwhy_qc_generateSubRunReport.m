function fmrwhy_qc_generateSubRunReport(bids_dir, sub, task, run, options)


% ------------------------------------
% STEP 1 -- Load default filenames etc
% ------------------------------------

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% TODO: have to loop through runs for this one if we want to generate report for all runs for sub:
% Update workflow params with subject functional derivative filenames
% Choose arbitrary echo for now, since this is not needed for current qc report
echo = '2';
options = fmrwhy_defaults_subFunc(bids_dir, sub, '', task, run, echo, options);

% ---------------------------------------------------------------
% STEP 2 -- Create structure with content to write to HTML report
% ---------------------------------------------------------------

report.param_asset_css = fullfile(options.fmrwhy_dir, 'assets/fmrwhy_reports.css');
report.param_asset_js = fullfile(options.fmrwhy_dir, 'assets/fmrwhy_reports.js');
report.param_avatar = fullfile(options.fmrwhy_dir, 'img/logo_jsheunis_3.jpeg');
report.param_sub = sub;
report.param_datetime = '2020 blabla';
report.param_brain_mask = fullfile(options.anat_dir_qc, ['sub-' sub '_brain_mask_montage.png']); % /Volumes/Stephan_WD/NEUFEPME_data_BIDS/derivatives/fmrwhy-qc/sub-001/anat/sub-001_brain_mask_montage.png
report.param_gm_mask = fullfile(options.anat_dir_qc, ['sub-' sub '_GM_mask_montage.png']);
report.param_wm_mask = fullfile(options.anat_dir_qc, ['sub-' sub '_WM_mask_montage.png']);
report.param_csf_mask = fullfile(options.anat_dir_qc, ['sub-' sub '_CSF_mask_montage.png']);
report.param_mean_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_mean.png']);
report.param_tsnr_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_tsnr.png']);
report.param_var_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_var.png']);
report.param_std_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_std.png']);
report.param_rograyplot_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-RO_grayplot.png']);
report.param_gsograyplot_img = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-GSO_grayplot.png']);
report.param_physioqc_img = fullfile(options.func_dir_qc, 'PhysIO_output_level2_04.jpg');


% ---------------------------
% STEP 3 -- Set up HTML files
% ---------------------------

% Set up datetime strings for filename and content
dt = datetime('now');
[Y,MO,D,H,MI,S] = datevec(dt);
dt_str = [num2str(Y) num2str(MO) num2str(D) num2str(H) num2str(MI) num2str(round(S))];
t = datestr(dt);

% Set up filenames
% TODO: store this in defaults somewhere, also need to be able to infer fMRwhy directory automatically / or run init function
html_template_fn = '/Users/jheunis/Documents/MATLAB/fMRwhy/fmrwhy_template_subQCreport.htm';
html_new_fn = fullfile(options.sub_dir_qc, ['sub-' sub '_desc-QC' dt_str '_report.html']);
html_tmp_fn = 'temp.html';


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
        fprintf(fid_tmp,'%s', newString);
    end
end
fclose(fid_template);
fclose(fid_tmp);

[status, msg, msgID] = movefile(html_tmp_fn, html_new_fn);
delete(html_tmp_fn)









%% Get sample HTML file
%fid1 = fopen(html_template_fn);
%lines = textscan(fid,'%s','delimiter','\n');
%fclose(fid);
%lines = lines{1};
%%replace tags with actual values
%relevant =  find(~cellfun(@isempty,strfind(s{1},'E_FROM')));
%for i = 1:length(relevant)
%    lines{relevant(i)-14} = '<modified line>';
%end
%% Create subject report
%fid = fopen('out.txt','w');
%for i = 1:length(lines)
%    fprintf(fid,'%s\n',lines{i});
%end
%fclose(fid)


























