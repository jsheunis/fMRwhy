function stats = fmrwhy_qc_createStatsOutput(bids_dir, sub, ses, task, run, echo, opts)
%function stats = fmrwhy_qc_calculateStats(bids_dir, sub, ses, task, run)
% Function to calculate multiple statistical measures from single run of (masked) fMRI timeseries data

% Setup fmrwhy bids directories on workflow level
fmrwhy_defaults_setupDerivDirs(bids_dir);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
fmrwhy_defaults_setupSubDirs(bids_dir, sub);

% Update workflow params with subject anatomical derivative filenames
opts = fmrwhy_defaults_subAnat(bids_dir, sub, opts);

% Update workflow params with subject functional derivative filenames
opts = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, opts);


% Stats images filenames
mean_fn = fullfile(opts.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_mean.nii']);
std_fn = fullfile(opts.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_std.nii']);
tsnr_fn = fullfile(opts.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_tsnr.nii']);
var_fn = fullfile(opts.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_var.nii']);
% Stats timeseries filename
stats_timeseries_fn = fullfile(opts.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_desc-stats_timeseries.tsv']);
% Stats summary filename
stats_summary_fn = fullfile(opts.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_desc-stats_summary.tsv']);

% Create stats images and files if they don't exist yet
stats_fn = {mean_fn, std_fn, tsnr_fn, var_fn, stats_timeseries_fn, stats_summary_fn};
run_stats = 0;
for i = 1:numel(stats_fn)
    if ~exist(stats_fn{i}, 'file')
        disp(['Stats output does not exist yet: ' stats_fn{i}]);
        run_stats = 1;
    end
end
if run_stats
    disp('Computing and saving statistical images and files for fMRI timeseries data')
    % Get stats for fMRI timeseries
    stats = fmrwhy_qc_calculateStats(bids_dir, sub, opts.rfunctional_fn); % TODO: decide which timeseries to use, processed or not
    % save images to file: fmrwhy_util_saveNifti(template_fn, img, new_fn, descrip, pinfo)
    fmrwhy_util_saveNifti(mean_fn, stats.data_3D_mean, opts.template_fn, 'fMRI timeseries mean', 0)
    fmrwhy_util_saveNifti(std_fn, stats.data_3D_stddev, opts.template_fn, 'fMRI timeseries standard deviation', 1)
    fmrwhy_util_saveNifti(tsnr_fn, stats.data_3D_tsnr, opts.template_fn, 'fMRI timeseries tSNR', 1)
    fmrwhy_util_saveNifti(var_fn, stats.data_3D_var, opts.template_fn, 'fMRI timeseries variance', 0)
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

    disp('Complete!')
    disp('---')
else
    disp('Statistical images and files already exist for fMRI timeseries data. Loading them now.')
    disp('---')
    stats.data_3D_mean = spm_read_vols(spm_vol(mean_fn));
    stats.data_3D_stddev = spm_read_vols(spm_vol(std_fn));
    stats.data_3D_tsnr = spm_read_vols(spm_vol(tsnr_fn));
    stats.data_3D_var = spm_read_vols(spm_vol(var_fn));
end

% Get screen size for plotting
scr_size = get(0,'ScreenSize');
dist = scr_size(4);
if scr_size(3) < dist
    dist = scr_size(3);
end

% Generate montage images, if they don't exist
stats_image_fns = {mean_fn, std_fn, var_fn, tsnr_fn};
stats_image_txt = {'mean', 'std', 'var', 'tsnr'};
stats_image_fieldname = {'data_3D_mean', 'data_3D_stddev', 'data_3D_var', 'data_3D_tsnr'};
stats_image_colormaps = {'gray', 'parula', 'parula', 'hot'};
for i = 1:numel(stats_image_fns)
    stats_montage_fns{i} = fullfile(opts.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_space-individual_' stats_image_txt{i} '.png']);
    %if ~exist(stats_montage_fns{i}, 'file')
        montage = fmrwhy_util_createMontage(stats.(stats_image_fieldname{i}), 9, 1, stats_image_txt{i}, stats_image_colormaps{i}, 'off', 'max');
        ax = montage.ax;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset; 
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
        set(ax,'xtick',[])
        set(ax,'xticklabel',[])
        set(ax,'ytick',[])
        set(ax,'yticklabel',[])
        set(ax,'ztick',[])
        set(ax,'zticklabel',[])
        print(montage.f,stats_montage_fns{i},'-dpng', '-r0')
    %end
end



