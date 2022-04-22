function options = fmrwhy_bids_setupStatsSubDirs(bids_dir, sub, options)

    if isempty(options) || nargin < 3
        options = struct;
    end

    % Derivatives directories
    options.deriv_dir = fullfile(bids_dir, 'derivatives');
    options.stats_dir = fullfile(options.deriv_dir, 'fmrwhy-stats');

    % Subject directories
    options.sub_dir_stats = fullfile(options.stats_dir, ['sub-' sub]);

    % Create, if necessary
    if ~exist(options.sub_dir_stats, 'dir')
        % Create new stats deriv dir
        mkdir(options.sub_dir_stats);
    end