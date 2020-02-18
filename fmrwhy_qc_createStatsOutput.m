function stats = fmrwhy_qc_createStatsOutput(bids_dir, sub, ses, task, run)
%function stats = fmrwhy_qc_calculateStats(bids_dir, sub, ses, task, run)
% Function to calculate multiple statistical measures from single run of (masked) fMRI timeseries data

% Load/create required defaults
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
template_task = 'rest';
template_run = '1';
template_echo = '2';

% Directory and content setup
deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
qc_dir = fullfile(deriv_dir, 'fmrwhy-qc');
sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);
sub_dir_qc = fullfile(qc_dir, ['sub-' sub]);
func_dir_qc = fullfile(sub_dir_qc, 'func');
if ~exist(func_dir_qc, 'dir')
    mkdir(func_dir_qc)
end

% Grab anatomical image
anatomical_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_T1w.nii']);

% Grab functional timeseries filename,
functional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_bold.nii']);

% Grab template filename
template_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_space-individual_bold.nii']);

% Stats images filenames
mean_fn = fullfile(func_dir_qc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_space-individual_mean.nii']);
std_fn = fullfile(func_dir_qc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_space-individual_std.nii']);
tsnr_fn = fullfile(func_dir_qc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_space-individual_tsnr.nii']);
var_fn = fullfile(func_dir_qc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_space-individual_var.nii']);
% Stats timeseries filename
stats_timeseries_fn = fullfile(func_dir_qc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_desc-stats_timeseries.tsv']);
% Stats summary filename
stats_summary_fn = fullfile(func_dir_qc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_desc-stats_summary.tsv']);

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
    stats = fmrwhy_qc_calculateStats(bids_dir, sub, functional_fn);
    % save images to file: fmrwhy_util_saveNifti(template_fn, img, new_fn, descrip, pinfo)
    fmrwhy_util_saveNifti(mean_fn, stats.data_3D_mean, template_fn, 'fMRI timeseries mean', 1)
    fmrwhy_util_saveNifti(std_fn, stats.data_3D_stddev, template_fn, 'fMRI timeseries standard deviation', 1)
    fmrwhy_util_saveNifti(tsnr_fn, stats.data_3D_tsnr, template_fn, 'fMRI timeseries tSNR', 1)
    fmrwhy_util_saveNifti(var_fn, stats.data_3D_var, template_fn, 'fMRI timeseries variance', 1)
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
    disp('Statistical images and files already exist for fMRI timeseries data')
    disp('---')
end


% Generate montage images