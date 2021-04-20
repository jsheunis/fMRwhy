function options = fmrwhy_settings_load(settings_fn)

    % --------------------------------------------------------------------------

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
    options = fmrwhy_bids_setupQcDerivDirs(options.bids_dir, options);

    % Validate settings
    options = fmrwhy_settings_validate(options);
