function fmrwhy_workflow_1stlevel_singleTask(bids_dir, sub, ses, task, options)
    % A custom workflow that runs 1st level analysis for all runs of a single task of
    % a single subject based on specified parameters.
    % This workflow assumes that an fMRwhy preprocessing workflow has been completed successfully
    % and that the preproc and qc derivatives are in place.
    %
    % :param bids_dir: Location of BIDS dataset
    % :type bids_dir: character array - path to directory
    % :param sub: BIDS subject code
    % :type sub: character array
    % :param ses: BIDS session code
    % :type ses: character array
    % :param tasks: BIDS task code
    % :type tasks: character array
    % :param options: the structure with parameter settings pertaining to the analysis and dataset
    % :type options: structure
    %
    % :Outputs: A directory with ...
    %
    % :Example:
    %
    % >>> 
    % >>> 

    % STEPS:
    % 0 - setup
    % 1 - loop through all runs for specified sub+ses+task, and
    % 1.1: grab preprocessed functional run +
    % 1.2: create multiple regressor file +
    % 1.3: set up statistical design parameters +


    % -------------------
    % SETUP (DO NOT EDIT)
    % -------------------
    % % SETUP STEP A -- Set fMRwhy defaults, check dependencies, set Matlab path, etc
    % options = fmrwhy_defaults();
    % options = fmrwhy_util_checkDependencies(options);
    % % SETUP STEP B -- Load settings, directories, filenames and parameters
    % % Run settings file ==> populates study/data-specific fields in the options structure, including BIDS variables
    % run(settings_fn);
    % % Setup fmrwhy derivative directories on workflow level
    % options = fmrwhy_bids_setupQcDerivDirs(options.bids_dir, options);
    % % Validate settings
    % options = fmrwhy_settings_validate(options);

    % ------------------
    % 1ST LEVEL PIPELINE
    % ------------------

    % -------
    % STEP 0: get runs for current task
    % -------
    runs = bids.query(options.bids_dataset, 'runs', 'sub', sub, 'ses', ses, 'task', task);
    
    % -------
    % STEP 1: set up directories for outputs, navigate to current task directory
    % -------
    options = fmrwhy_bids_setupStatsTaskDirs(options.bids_dir, sub, ses, task, options);
    task_dir_stats = options.task_dir_stats;
    cd(task_dir_stats);

    % -------
    % STEP 2: Loop through runs and:
    % -------
    glm_settings = options.firstlevel.glm_regressors;
    timing_params_task = options.firstlevel.(task).timing_params;
    sess_params_task = options.firstlevel.(task).sess_params;
    cond_params_task = options.firstlevel.(task).cond_calcs;
    confounds_fn = {};
    regressors_fn = {};
    sess_params = {};
    functional_fn = {};
    regressors_names = {};
    for r = 1:numel(runs)
        rn = runs{r};
        % Update options structure for current run
        options = fmrwhy_bids_getFuncDerivs(options.bids_dir, sub, task, options, 'ses', ses, 'task', task, 'run', rn);

        % -------
        % STEP 2.1: grab preprocessed functional run
        % -------
        [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', sub, 'ses', ses, 'task', task, 'run', rn, 'desc', 'swrapreproc', 'ext', '_bold.nii');
        functional_fn{r} = fullfile(options.preproc_dir, filepath, filename);

        % -------
        % STEP 2.2: create multiple regressor file
        % -------
        confounds_fn{r} = options.confounds_fn;
        regressors_fn{r} = fullfile(task_dir_stats, ['sub-' sub, '_ses-' ses, '_task-' task '_run-' rn '_desc-GLM_regressors.txt']);
        [regressors_mat, regressors_names{r}] = fmrwhy_1stlevel_createRegressors(confounds_fn{r}, regressors_fn{r}, glm_settings);

        % -------
        % STEP 2.3: set up statistical design parameters
        % -------
        [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', sub, 'ses', ses, 'task', task, 'run', rn, 'ext', '_events.tsv');
        event_filename = fullfile(options.preproc_dir, filepath , filename);

        sess_params{r} = sess_params_task;
        sess_params{r}.cond = fmrwhy_1stlevel_filterCalcEvents(event_filename, sess_params{r}.cond_names, cond_params_task);
        % [cond, trials] = fmrwhy_util_1stlevelBIDStoConditions(event_filename, sess_params{r}.cond_names);
        % sess_params{r}.cond = cond;
    end

    % -------
    % STEP 3: Create the 1st level model
    % -------
    fmrwhy_batch_specify1stlevel(task_dir_stats, functional_fn, regressors_fn, sess_params, timing_params_task);
    load([task_dir_stats filesep 'SPM.mat']);

    % -------
    % STEP 4: Estimate the 1st level model
    % -------
    fmrwhy_batch_estimate1stlevel(task_dir_stats);

    % -------
    % STEP 5: Review the model and output diagnostic figures (review is done once per output figure type)
    % -------
    display_options = options.firstlevel.review.display_options;
    output_filenames = options.firstlevel.review.output_filenames;
    for i = 1:numel(display_options)
        review_params = options.firstlevel.review.review_params
        review_params.display = display_options{i};
        fmrwhy_batch_review1stlevel(task_dir_stats, review_params);
        new_path = fmrwhy_1stlevel_renameReviewOutputs(task_dir_stats, output_filenames{i})
    end

    % -------
    % STEP 6: Setup and run contrasts
    % -------
    assignin('base','SPM',SPM)
    assignin('base','options',options)
    contrast_params = options.firstlevel.(task).contrast_params;
    consess = fmrwhy_1stlevel_extractContrasts(SPM, contrast_params, runs);
    fmrwhy_batch_contrast1stlevel(task_dir_stats, consess);

    % -------
    % STEP 7: Thresholding, results report, output figures
    % -------
    conspec_params = options.firstlevel.conspec_params;
    conspec = struct;
    for k = 1:numel(consess)
        conspec(k).titlestr = consess{k}.tcon.name;
        conspec(k).contrasts = k;
        conspec(k).threshdesc = conspec_params.threshdesc;
        conspec(k).thresh = conspec_params.thresh;
        conspec(k).extent = conspec_params.extent;
        conspec(k).conjunction = conspec_params.conjunction;
        conspec(k).mask = conspec_params.mask;
    end
    fmrwhy_batch_threshold1stlevel(task_dir_stats, conspec);
    % Rename files
    dt = datetime;
    y = num2str(year(dt));
    m = month(dt, 'shortname');
    m = m{1};
    d = sprintf('%02d', day(dt));
    for k = 1:numel(consess)
        src = fullfile(task_dir_stats, ['spm_' y m d '_' sprintf('%03d', k) '.jpg']);
        dest = fullfile(task_dir_stats, ['statsresults_clusters_' sprintf('%03d', k) '.jpg']);
        movefile(src, dest);
    end

    % TODO: how to get xSPM without the code below, which for some reason prompts SPM dialog box to open up and ask for file selection
    % [SPM, xSPM] = spm_getSPM(fullfile(run_dir_stats, 'SPM.mat'));
    % assignin('base', 'SPM', SPM)
    % if exist('xSPM','var')
    %    disp('xSPM exists. Saving.')
    %    save('xSPM.mat', 'xSPM')
    % else
    %    disp('xSPM does not exist.')
    % end

    % -------
    % STEP 9: Tmap montages
    % -------
    % TODO: add option to detect standard or subject space, and automatically grab the correct background file
    canon_fn = fullfile(options.spm_dir, 'canonical/single_subj_T1.nii');
    canon_stat_fn = fullfile(task_dir_stats, 'temp_single_subj_T1.nii');
    copyfile(canon_fn, canon_stat_fn);
    template_fn = fullfile(task_dir_stats, 'mask.nii');
    fmrwhy_batch_realignResl({template_fn, canon_stat_fn});
    background_fn = fullfile(task_dir_stats, 'rrtemp_single_subj_T1.nii');
    % background_fn = options.rcoregest_anatomical_fn;
    [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(background_fn);
    background_img = p.nii.img;

    for k = 1:numel(consess)
        tmap_fn = fullfile(task_dir_stats, ['spmT_' sprintf('%04d', k) '.nii']);
        tmap_clusters_fn = fullfile(task_dir_stats, ['spmT_' sprintf('%04d', k) '_binary_clusters.nii']);
        [ptmap, ~, ~, ~] = fmrwhy_util_readOrientNifti(tmap_fn);
        [ptmapc, ~, ~, ~] = fmrwhy_util_readOrientNifti(tmap_clusters_fn);
        stats_img = fmrwhy_util_maskImage(double(ptmap.nii.img), double(ptmapc.nii.img));
        str = consess{k}.tcon.name;
        saveAs_fn = fullfile(task_dir_stats, ['sub-' sub '_ses-' ses '_task-' task '_desc-' str '_threshtmap.png']);
        overlaymontage = fmrwhy_util_createStatsOverlayMontage(p.nii.img, stats_img, [], 12, 1, '', 'gray', 'off', 'max', [], 'hot', [], true, saveAs_fn);
    end

    clear SPM