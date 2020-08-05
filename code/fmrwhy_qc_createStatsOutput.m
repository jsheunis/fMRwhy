function stats = fmrwhy_qc_createStatsOutput(bids_dir, sub, ses, task, run, echo, options)
%function stats = fmrwhy_qc_calculateStats(bids_dir, sub, ses, task, run)
% Function to calculate multiple statistical measures from single run of (masked) fMRI timeseries data

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


% Stats images filenames
mean_fn = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_mean.nii']);
std_fn = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_std.nii']);
tsnr_fn = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_tsnr.nii']);
var_fn = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_var.nii']);

% Stats timeseries filename
stats_timeseries_fn = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_desc-stats_timeseries.tsv']);
% Stats summary filename
stats_summary_fn = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_desc-stats_summary.tsv']);

% Create stats images and files if they don't exist yet
stats_fn = {mean_fn, std_fn, tsnr_fn, var_fn, stats_timeseries_fn, stats_summary_fn};
run_stats = 0;
stats = [];
for i = 1:numel(stats_fn)
    if ~exist(stats_fn{i}, 'file')
        disp(['Stats output does not exist yet: ' stats_fn{i}]);
        run_stats = 1;
    end
end
if run_stats || options.qc_overwrite_statsoutput
    disp('Computing and saving statistical images and files for fMRI timeseries data')
    % Get stats for fMRI timeseries
    stats = fmrwhy_qc_calculateStats(bids_dir, sub, options.rafunctional_fn); % TODO: decide which timeseries to use, processed or not, Aug 4 2020 note: changed this from rfunctional to rafunctional when redoing QC output. Reasoning = to be consistent with multi-echo use.
    % save images to file: fmrwhy_util_saveNifti(template_fn, img, new_fn)
    no_scaling = 1;
    fmrwhy_util_saveNifti(mean_fn, double(stats.data_3D_mean), options.template_fn, no_scaling);
    fmrwhy_util_saveNifti(std_fn, double(stats.data_3D_stddev), options.template_fn, no_scaling);
    fmrwhy_util_saveNifti(tsnr_fn, double(stats.data_3D_tsnr), options.template_fn, no_scaling);
    fmrwhy_util_saveNifti(var_fn, double(stats.data_3D_var), options.template_fn, no_scaling);
    % Write stats timeseries data to tsv file
    [d, f, e] = fileparts(stats_timeseries_fn);
    temp_txt_fn = fullfile(d, [f '.txt']);
    col_names = {'zscore'};
    data = stats.data_2D_zstat_ts;
    data_table = array2table(data,'VariableNames', col_names);
    writetable(data_table, temp_txt_fn, 'Delimiter','\t');
    [status, msg, msgID] = movefile(temp_txt_fn, stats_timeseries_fn);
    % Write stats summary data to tsv file
    [d, f, e] = fileparts(stats_summary_fn);
    temp_txt_fn = fullfile(d, [f '.txt']);
    col_names = {'tSNR_mean_GM', 'tSNR_mean_WM', 'tSNR_mean_CSF', 'tSNR_mean_brain', 'zscore_mean', 'gcor'};
    data = [stats.tSNR_mean_GM stats.tSNR_mean_WM stats.tSNR_mean_CSF stats.tSNR_mean_brain stats.Zstat_mean stats.GCOR];
    data_table = array2table(data,'VariableNames', col_names);
    writetable(data_table, temp_txt_fn, 'Delimiter','\t');
    [status, msg, msgID] = movefile(temp_txt_fn, stats_summary_fn);

    % Generate montage images, if they don't exist
    imgs = struct;
    imgs_masked = struct;
    stats_image_fns = {mean_fn, std_fn, var_fn, tsnr_fn};
    stats_image_txt = {'mean', 'std', 'var', 'tsnr'};
    stats_image_fieldname = {'data_3D_mean', 'data_3D_stddev', 'data_3D_var', 'data_3D_tsnr'};
    stats_image_colormaps = {'gray', 'parula', 'parula', 'hot'};
    stats_image_cxs = {[], [], [], [0 250]};
    stats_image_clrbr = {false, false, false, true};

    % Plotting settings
    rgb_ongray = [255, 115, 236];
    rgb_onhot = [148, 239, 255];
    rgb_onparula = [255, 115, 236];
    use_rgb = {rgb_ongray, rgb_onparula, rgb_onparula, rgb_onhot};

    % Mask settings
    masks_oriented = fmrwhy_util_loadOrientMasks(bids_dir, sub);
    mask_img_oriented = masks_oriented.brain_mask_3D;

    for i = 1:numel(stats_image_fns)
        stats_montage_fns{i} = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_' stats_image_txt{i} '.png']);
        stats_montage_masked_fns{i} = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_desc-brainMasked_' stats_image_txt{i} '.png']);
        if ~exist(stats_montage_fns{i}, 'file') || options.qc_overwrite_statsoutput

            [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(stats_image_fns{i});

            imgs.(stats_image_fieldname{i}) = double(p.nii.img);
            imgs_masked.(stats_image_fieldname{i}) = fmrwhy_util_maskImage(double(p.nii.img), mask_img_oriented);

            montage = fmrwhy_util_createStatsOverlayMontage(imgs.(stats_image_fieldname{i}), [], [], 9, 1, stats_image_txt{i}, stats_image_colormaps{i}, 'off', 'max', stats_image_cxs{i}, [], use_rgb{i}, stats_image_clrbr{i}, stats_montage_fns{i});
            montage_masked = fmrwhy_util_createStatsOverlayMontage(imgs_masked.(stats_image_fieldname{i}), [], [], 9, 1, stats_image_txt{i}, stats_image_colormaps{i}, 'off', 'max', stats_image_cxs{i}, [], use_rgb{i}, stats_image_clrbr{i}, stats_montage_masked_fns{i});
            % fmrwhy_util_createStatsOverlayMontage(template_img, stats_img, roi_img, columns, rotate, str, clrmp, visibility, shape, cxs, stats_clrmp, roi_rgbcolors, clrbar, saveAs_fn)
        end
    end

    disp('Complete!')
    disp('---')
else
    disp('Statistical images and files already exist for fMRI timeseries data.')
    disp('---')
end



