function fmrwhy_defaults_setupSubDirs(bids_dir, sub)

% Derivatives directories
deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
qc_dir = fullfile(deriv_dir, 'fmrwhy-qc');
stats_dir = fullfile(deriv_dir, 'fmrwhy-stats');

% Subject directories
sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);
sub_dir_qc = fullfile(qc_dir, ['sub-' sub]);
sub_dir_stats = fullfile(stats_dir, ['sub-' sub]);
sub_dir_BIDS = fullfile(bids_dir, ['sub-' sub]);

% Create and copy content, if necessaru
if ~exist(sub_dir_preproc, 'dir')
    mkdir(sub_dir_preproc)
    copyfile(sub_dir_BIDS, sub_dir_preproc)
end
if ~exist(sub_dir_qc, 'dir')
    mkdir(sub_dir_qc)
end
if ~exist(sub_dir_stats, 'dir')
    mkdir(sub_dir_stats)
end