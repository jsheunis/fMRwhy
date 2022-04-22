function options = fmrwhy_bids_setupStatsTaskDirs(bids_dir, sub, ses, task, options)

    if isempty(options) || nargin < 5
        options = struct;
    end

    options = fmrwhy_bids_setupStatsDerivDirs(bids_dir, options);
    options = fmrwhy_bids_setupStatsSubDirs(bids_dir, sub, options);

    if isempty(ses)
        task_dir_stats = fullfile(options.sub_dir_stats, task);
    else
        task_dir_stats = fullfile(options.sub_dir_stats, ses, task);
    end
    
    if ~exist(task_dir_stats, 'dir')
        mkdir(task_dir_stats);
    end

    options.task_dir_stats = task_dir_stats;