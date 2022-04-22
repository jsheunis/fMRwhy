function new_path = fmrwhy_1stlevel_renameReviewOutputs(stats_dir, new_name)
    % Rename the figures output by ``fmrwhy_batch_review1stlevel(task_dir_stats, review_params)``
    %
    % :param stats_dir: The directory within which the ``fmrwhy_batch_review1stlevel`` command was run
    % :type stats_dir: string or character array
    % :param new_name: The (partial) new name for the file
    % :type new_name: string or character array
    % :returns: ``new_path`` - path(s) to renamed file(s)

    % Rename figure outputs
    dt = datetime;
    y = num2str(year(dt));
    m = month(dt, 'shortname');
    m = m{1};
    d = sprintf('%02d', day(dt));
    fnames = dir([stats_dir filesep 'spm_' y m d '_*.jpg']);
    if isempty(fnames)
        disp('Warning: stats review file outputs not found... please inspect');
    else
        new_path = {};
        for j = 1:numel(fnames)
            src = fullfile(stats_dir, ['spm_' y m d '_' sprintf('%03d', j) '.jpg']);
            dest = fullfile(stats_dir, ['statsreview_' new_name '_' sprintf('%03d', j) '.jpg']);
            movefile(src, dest);
            new_path = [new_path, {dest}];
        end
    end

    
    