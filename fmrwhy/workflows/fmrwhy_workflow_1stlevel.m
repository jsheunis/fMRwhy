function options = fmrwhy_workflow_1stlevel(settings_fn)
    % A custom workflow that runs 1st level analysis for all subjects in a BIDS directory.
    % This workflow assumes that an fMRwhy preprocessing workflow has been completed successfully
    % and that the preproc and qc derivatives are in place.
    %
    % :param settings_fn: User preference settings for the specific analysis and dataset, see
    % :type settings_fn: character array - path to M-file
    % :returns: ``options`` - the updated structure with parameter settings pertaining to the analysis and dataset
    %
    % :Outputs: A directory with ...
    %
    % :Example:
    %
    % >>> settings_fn = '/fMRwhy/fmrwhy/settings/fmrwhy_settings_template.m';
    % >>> options = fmrwhy_workflow_1stlevel(settings_fn)
    % 
    % STEPS:
    % 0 - setup
    % 1 - loop through all subjects (from settings) and call fmrwhy_workflow_1stlevel_singleSub

    options = fmrwhy_defaults();

    % -------
    % SETUP STEP A -- Check dependencies, Matlab path, etc
    % -------
    options = fmrwhy_util_checkDependencies(options);

    % -------
    % SETUP STEP B -- Load settings, defaults, filenames and parameters
    % -------
    % Run settings file ==> populates study/data-specific fields in the options structure, including BIDS variables
    run(settings_fn);

    % Setup fmrwhy derivative directories on workflow level
    options = fmrwhy_bids_setupStatsDerivDirs(options.bids_dir, options);

    % Validate settings
    options = fmrwhy_settings_validate(options);

    % Load the subjects, sessions, tasks, runs
    % subs = options.subjects_output;
    subs = {'mcc000701'};
    % sessions = options.sessions;
    sessions = {'w03lab'};
    % tasks = options.tasks;
    tasks = {'snat'};
    runs = options.runs;
    % TODO: check for and handle empty cell arrays here
    % TODO: decide on a strategy for handling varying sessions/tasks/runs per subject

    % -------
    % 1st level Pipeline -- for each subject
    % -------
    disp('------------------------');
    disp(['START 1st LEVEL ANALYSIS for dataset at: ' options.bids_dir]);
    disp('------------------------');
    for s = 1:numel(subs)
        sub = subs{s};

        % -------
        % STEP 1 -- 1st level analysis
        % -------
        fmrwhy_workflow_1stlevel_singleSub(sub, sessions, tasks, runs, settings_fn, options.spm_dir)
    end
