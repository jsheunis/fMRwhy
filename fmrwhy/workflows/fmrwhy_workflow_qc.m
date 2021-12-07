function options = fmrwhy_workflow_qc(settings_fn)
    % A custom workflow that does anatomical and functional data quality control for all subjects in a BIDS directory.
    % Steps include anatomical-to-functional processing, basic functional time series preprocessing,
    % generating a range of QC metrics and images, and compiling an HTML-report per subject.
    %
    % :param settings_fn: User preference settings for the specific analysis and dataset, see
    % :type settings_fn: M-file
    % :returns: ``options`` - the updated structure with parameter settings pertaining to the analysis and dataset
    %
    % :Outputs: A directory with all assets necessary for rendering the HTML report, located in ``derivatives/fmrwhy-qc/sub-XXX/report_[yyyymmddhhmmss]``
    %
    % :Example:
    %
    % >>> settings_fn = '/fMRwhy/fmrwhy/settings/fmrwhy_settings_template.m';
    % >>> options = fmrwhy_workflow_qc(settings_fn)
    %
    % .. seealso::
    %   - A thorough description of the :ref:`quality_reporting` pipeline
    %   - A sample report can be `viewed here`_.
    %
    % .. _viewed here: https://jsheunis.github.io/fmrwhy_sample_QCreport.html

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

    % Load the subjects, sessions, tasks, runs
    subs = options.subjects_output;
    sessions = options.sessions;
    tasks = options.tasks;
    runs = options.runs;
    % TODO: check for and handle empty cell arrays here
    % TODO: decide on a strategy for handling varying sessions/tasks/runs per subject

    % -------
    % QC Pipeline -- for each subject
    % -------

    for s = 1:numel(subs)
        sub = subs{s};

        % -------
        % STEP 1 -- Preprocessing and QC steps
        % -------
        fmrwhy_workflow_qc_singleSub(sub, sessions, tasks, runs, settings_fn, options.spm_dir)

        % -------
        % STEP 2 -- QC report: fmrwhy_workflow_qcSubReport.m
        % -------
        fmrwhy_workflow_qcSubReport(sub, options);

    end
