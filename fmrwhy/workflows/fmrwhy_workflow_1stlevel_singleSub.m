function fmrwhy_workflow_1stlevel_singleSub(sub, sessions, tasks, runs, settings_fn, spm_dir)
    % A custom workflow that runs 1st level analysis for a single subject based on specified parameters
    % This workflow assumes that an fMRwhy preprocessing workflow has been completed successfully
    % and that the preproc and qc derivatives are in place.
    %
    % :param sub: BIDS subject code
    % :type sub: character array
    % :param sessions: BIDS sessions
    % :type sessions: cell array of character arrays
    % :param tasks: BIDS tasks
    % :type tasks: cell array of character arrays
    % :param runs: BIDS runs
    % :type runs: cell array of character arrays
    % :param settings_fn: User preference settings for the specific analysis and dataset, see
    % :type settings_fn: character array - path to M-file
    % :param spm_dir: Directory where SPM12 is installed
    % :type spm_dir: character array - path to SPM12 directory
    % 
    % :returns: ``options`` - the updated structure with parameter settings pertaining to the analysis and dataset
    %
    % :Outputs: A directory with ...
    %
    % :Example:
    %
    % >>> 
    % >>> 

    % STEPS:
    % 0 - setup
    % 1 - loop through all sessions (argument), and
    % 2 - loop through all tasks (argument), and call fmrwhy_workflow_1stlevel_singleTask


    % -------------------
    % SETUP (DO NOT EDIT)
    % -------------------
    % SETUP STEP A -- Set fMRwhy defaults, check dependencies, set Matlab path, etc
    options = fmrwhy_defaults();
    options = fmrwhy_util_checkDependencies(options);
    % SETUP STEP B -- Load settings, directories, filenames and parameters
    % Run settings file ==> populates study/data-specific fields in the options structure, including BIDS variables
    run(settings_fn);
    % Setup fmrwhy derivative directories on workflow level
    options = fmrwhy_bids_setupStatsDerivDirs(options.bids_dir, options);
    options = fmrwhy_bids_getStatsDerivDirs(options.bids_dir, options);
    % Validate settings
    options = fmrwhy_settings_validate(options);

    % ------------------
    % 1ST LEVEL PIPELINE
    % ------------------

    % The following is executed per session, and per task.
    % For runs, the steps are either executed on the task level (i.e. all runs in one step),
    % or on the individual run level, depending on the specific step.
    % TODO: this currently assumes that there are sessions, need to fully address the possibility of no sessions.

    % Setup fmrwhy derivatives directories on subject level
    options = fmrwhy_bids_setupStatsSubDirs(options.bids_dir, sub, options);

    for ss = 1:numel(sessions)
        ses = sessions{ss};

        % Skip session loop if current session does not exist
        ses_dir = fullfile(options.preproc_dir, ['sub-' sub], ['ses-' ses]);
        if ~exist(ses_dir, 'dir')
            disp('---');
            disp(['No session ses-', ses, ' for sub-', sub]);
            disp('---');
            continue;
        end

        % Skip session loop if there is no functional data for the current session
        func_dir = fullfile(options.preproc_dir, ['sub-' sub], ['ses-' ses], 'func');
        if ~exist(func_dir, 'dir')
            disp('---');
            disp(['No functional data for sub-', sub, '_ses-', ses]);
            disp('---');
            continue;
        end

        % Loop through tasks
        for t = 1:numel(tasks)
            task = tasks{t};

            % Check if files with current task are available; if not, skip iteration
            fnames = dir([func_dir filesep '*' task '*' ]);
            if isempty(fnames)
                disp('---');
                disp(['No functional files for sub-', sub, '_ses-', ses, '_task-', task]);
                disp('---');
                continue;
            end

            disp('---------------------');
            disp(['SUB:  ', sub]);
            disp(['SES:  ', ses]);
            disp(['TASK: ', task]);
            disp('---------------------');

            fmrwhy_workflow_1stlevel_singleTask(options.bids_dir, sub, ses, task, options)
        end
    end




% % Loop through runs
% for r = 1:numel(runs)
%     rn = runs{r};

%     % Check if files with current task and run are available; if not, skip iteration
%     f_run_names = dir([func_dir filesep '*' task '*' 'run-' run '*']);
%     if isempty(f_run_names)
%         disp('---');
%         disp(['No files for sub-', sub, '_ses-', ses, '_task-', task, '_run-', run]);
%         disp('---');
%         continue;
%     end

%     disp('---------------------');
%     disp(['START 1st LEVEL ANALYSIS: sub-', sub, '_ses-', ses, '_task-', task, '_run-', run]);
%     disp('---------------------');

%     % STEP 1: set up directories for outputs
%     options = fmrwhy_bids_setupStatsTaskRunDirs(options.bids_dir, sub, ses, task, run, options)
%     run_dir_stats = options.run_dir_stats;
%     % STEP 2: navigate to current session-task-run directory
%     cd(run_dir_stats);
%     % STEP 3: run-specific 1st level analysis
%     fmrwhy_workflow_1stlevel_singleRun(options.bids_dir, sub, ses, task, run, echo, options)
% end