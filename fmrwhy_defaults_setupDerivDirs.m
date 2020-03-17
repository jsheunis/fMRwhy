function fmrwhy_defaults_setupDerivDirs(bids_dir)

% Derivatives directory setup
deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
qc_dir = fullfile(deriv_dir, 'fmrwhy-qc');
stats_dir = fullfile(deriv_dir, 'fmrwhy-stats');

dirs = {preproc_dir, qc_dir, stats_dir};
for i = 1:numel(dirs)
    if ~exist(dirs{i}, 'dir')
        mkdir(dirs{i})
    end
end