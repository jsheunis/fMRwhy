function fmrwhy_workflow_1stlevel_singleRun(bids_dir, sub, ses, task, run, echo, options)

    % A custom workflow that runs 1st level analysis for a single subject based on specified parameters

    % STEPS:
    % 

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
    options = fmrwhy_bids_setupQcDerivDirs(options.bids_dir, options);
    % Validate settings
    options = fmrwhy_settings_validate(options);

    % -------------------------------------------
    % SUBJECTS, SESSIONS, TASKS, RUNS (DO NOT EDIT THIS, unless you specifically want to override defaults)
    % -------------------------------------------
    % Create a cell array with the required sessions to run
    % sessions = {'1', '2', '3'};
    % Create a cell array with the required tasks
    % tasks = {'rest', 'motor'};
    % Create a cell array with the required runs
    % runs = {'1', '2', '3'};

    % NOTE: for now, we're assuming that all subjects have all sessions, all sessions have all tasks, and all tasks has all runs.
    % During the workflow, we loop through all possible combinations of these,
    % and we do a check to see if the corresponding files actually exist before starting the run-specific pipeline.
    % If it doesn't exist, we write a notification to the command line and then continue with the next iteration.
    % There are smarter ways of building this logic, but this is the easiest for now. TODO.

    % ------------------
    % 1ST LEVEL PIPELINE
    % ------------------

    % Load/create required defaults
    % Setup fmrwhy BIDS-derivatuve directories on workflow level
    options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

    % Grab parameters from workflow settings file
    options = fmrwhy_settings_preprocQC(bids_dir, options);

    % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

    % set flag to specify echo information
    % If one of the standard echoes was passed, set flag to true.
    % If non-standard echo string was passed, set flag to false, assuming that this points to a string to identify a combined echo timeseries in the fmrwhy-multiecho derivatives folder
    echoes = {};
    is_standard_echo = false;
    for e = 1:options.Ne
        echoes{e} = num2str(e);
    end
    if sum(ismember(echoes, echo))
        is_standard_echo = true;
    end

    % Update workflow params with subject functional derivative filenames
    if is_standard_echo
        options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);
    else
        options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, options.template_echo, options);
    end

    options = fmrwhy_bids_getFuncDerivs(options.bids_dir, sub, task, options, 'ses', ses, 'task', task, 'run', rn);

    % -------
    % STEP 1: set up directories for outputs
    % -------
    options.sub_dir_stats = fullfile(options.stats_dir, ['sub-' sub]);
    run_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-' echo]);
    if ~exist(run_dir_stats, 'dir')
        mkdir(run_dir_stats);
    end
    cd(run_dir_stats);
    % -------
    % STEP 2: create design matrix regressors
    % -------
    % Load multiple confound regressors
    confounds_struct = tdfread(options.confounds_fn);
    confounds_mat = struct2array(confounds_struct);
    regressors_mat = [];
    regressors_names = {};

    fields = fieldnames(options.firstlevel.glm_regressors);

    for i = 1:numel(fields)
        key = fields{i};
        val = options.firstlevel.glm_regressors.(key);
        % If the value is true or larger than zero (for retroicor order), parse key and include relevant data in regressor matrix
        if val
            if strfind(key, 'trans_rot')
                if strcmp(key, 'trans_rot')
                    trans_rot_keys = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
                elseif strcmp(key, 'trans_rot_derivative1')
                    trans_rot_keys = {'trans_x_derivative1', 'trans_y_derivative1', 'trans_z_derivative1', 'rot_x_derivative1', 'rot_y_derivative1', 'rot_z_derivative1'};
                elseif strcmp(key, 'trans_rot_power2')
                    trans_rot_keys = {'trans_x_power2', 'trans_y_power2', 'trans_z_power2', 'rot_x_power2', 'rot_y_power2', 'rot_z_power2'};
                elseif strcmp(key, 'trans_rot_derivative1_power2')
                    trans_rot_keys = {'trans_x_derivative1_power2', 'trans_y_derivative1_power2', 'trans_z_derivative1_power2', 'rot_x_derivative1_power2', 'rot_y_derivative1_power2', 'rot_z_derivative1_power2'};
                else
                end
                for j = 1:numel(trans_rot_keys)
                    regressors_mat = [regressors_mat confounds_struct.(trans_rot_keys{j})];
                    regressors_names = [regressors_names {trans_rot_keys{j}}];
                end
            elseif strcmp(key, 'retroicor_c') || strcmp(key, 'retroicor_r') || strcmp(key, 'retroicor_cxr')
                for k = 1:val
                    new_key = [key num2str(k)];
                    regressors_mat = [regressors_mat confounds_struct.(new_key)];
                    regressors_names = [regressors_names {new_key}];
                end
            else
                regressors_mat = [regressors_mat confounds_struct.(key)];
                regressors_names = [regressors_names {key}];
            end
        end
    end

    regressors_fn = fullfile(run_dir_stats, ['sub-' sub '_task-' task '_run-' run '_echo-' options.template_echo '_desc-GLM_regressors.txt']);
    dlmwrite(regressors_fn, regressors_mat, 'delimiter', '\t', 'precision', '%1.7e');

    % -------
    % STEP 3: Set up statistical design parameters, based on task data
    % -------
    % Load task events file
    events_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_events.tsv']);
    events_struct = tdfread(events_fn);
    % Load session parameters (specifically, conditions, onsets and durations) from events file
    cond_names = options.firstlevel.(task).(['run' run]).sess_params.cond_names;
    [cond, trials] = fmrwhy_util_1stlevelBIDStoConditions(events_fn, cond_names);
    options.firstlevel.(task).(['run' run]).sess_params.cond = cond;
    sess_params = options.firstlevel.(task).(['run' run]).sess_params;
    % Select functional timeseries to use
    % If one of the standard echoes was passed as a parameter select appropriate timeseries, otherwise assume multiecho combined and select appropriate combined timeseries
    if is_standard_echo
        functional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srapreproc_bold.nii']);
    else
        options.me_dir = fullfile(options.deriv_dir, 'fmrwhy-multiecho');
        options.sub_dir_me = fullfile(options.me_dir, ['sub-' sub]);
        options.func_dir_me = fullfile(options.sub_dir_me, 'func');
        functional_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-s' echo '_bold.nii']);
    end

    % functional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-scombinedMEt2star_bold.nii'])

    % -------
    % STEP 4: Create the 1st level model
    % -------
    fmrwhy_batch_specify1stlevel(run_dir_stats, functional_fn, regressors_fn, sess_params);
    load([run_dir_stats filesep 'SPM.mat']);

    % -------
    % STEP 5: Estimate the 1st level model
    % -------
    fmrwhy_batch_estimate1stlevel(run_dir_stats);

    % -------
    % STEP 6: Review the model, output diagnostic figures
    % -------
    review_params.print = 'jpg';
    review_params_display = {'matrix', 'covariance', 'orth'};
    stats_review_outputs = {'matrix', 'covariance', 'orthogonality'};
    for i = 1:numel(review_params_display)
        review_params.display = review_params_display{i};
        fmrwhy_batch_review1stlevel(run_dir_stats, review_params);
        % Rename figure outputs
        dt = datetime;
        y = num2str(year(dt));
        m = month(dt, 'shortname');
        m = m{1};
        d = sprintf('%02d', day(dt));
        fnames = dir([run_dir_stats filesep 'spm_' y m d '_*.jpg']);
        if isempty(fnames)
            disp('Warning: stats review file outputs not found... please inspect');
            continue
        else
            %        bv1 = BVQXfile(fullfile(pth, fnames(1).name));
            for j = 1:numel(fnames)
                src = fullfile(run_dir_stats, ['spm_' y m d '_' sprintf('%03d', j) '.jpg']);
                dest = fullfile(run_dir_stats, ['statsreview_' stats_review_outputs{i} '_' sprintf('%03d', j) '.jpg']);
                movefile(src, dest);
            end
        end
    end

    % -------
    % STEP 7: Setup and run contrasts
    % -------
    [Ntt, Nregr] = size(SPM.xX.X);

    consess = options.firstlevel.(task).(['run' run]).contrast_params.consess;
    for j = 1:numel(consess)
        default_weights = consess{j}.tcon.weights;
        consess{j}.tcon.weights = zeros(1, Nregr);
        consess{j}.tcon.weights(1:length(default_weights)) = default_weights;
    end
    fmrwhy_batch_contrast1stlevel(run_dir_stats, consess);

    % -------
    % STEP 8: Thresholding, results report, output figures
    % -------
    conspec = struct;
    for k = 1:numel(consess)
        conspec(k).titlestr = consess{k}.tcon.name;
        conspec(k).contrasts = k;
        conspec(k).threshdesc = 'FWE';
        conspec(k).thresh = 0.0500;
        conspec(k).extent = 0;
        conspec(k).conjunction = 1;
        conspec(k).mask.none = 1;
    end
    fmrwhy_batch_threshold1stlevel(run_dir_stats, conspec);
    % Rename files
    for k = 1:numel(consess)
        src = fullfile(run_dir_stats, ['spm_' y m d '_' sprintf('%03d', k) '.jpg']);
        dest = fullfile(run_dir_stats, ['statsresults_clusters_' sprintf('%03d', k) '.jpg']);
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
    background_fn = options.rcoregest_anatomical_fn;
    [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(background_fn);
    background_img = p.nii.img;

    for k = 1:numel(consess)
        tmap_fn = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', k) '.nii']);
        tmap_clusters_fn = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', k) '_binary_clusters.nii']);
        [ptmap, ~, ~, ~] = fmrwhy_util_readOrientNifti(tmap_fn);
        [ptmapc, ~, ~, ~] = fmrwhy_util_readOrientNifti(tmap_clusters_fn);
        stats_img = fmrwhy_util_maskImage(double(ptmap.nii.img), double(ptmapc.nii.img));
        str = consess{k}.tcon.name;
        saveAs_fn = fullfile(run_dir_stats, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-' str '_threshtmap.png']);
        overlaymontage = fmrwhy_util_createStatsOverlayMontage(p.nii.img, stats_img, [], 9, 1, '', 'gray', 'off', 'max', [], 'hot', [], true, saveAs_fn);
    end

    % -------
    % STEP 10: Calculate anatomical-functional overlap ROI
    % -------
    % For now, just use the first ROI in the list of ROIs per task. TODO: change in future
    anat_roi_fn = {};
    anat_roi_fn{1} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-r' options.roi.(task).desc{1} '_roi.nii']);

    % For each anatomical ROI:
    for j = 1:numel(anat_roi_fn)
        % Grab binary image of ROI, not yet for plotting
        nii = nii_tool('load', anat_roi_fn{j});
        anat_roi_img = fmrwhy_util_createBinaryImg(nii.img, 0.1);
        % For each task contrast, i.e. for each tmap, grab binary tmap clusters, find overlap with anatomical ROI, save
        for k = 1:numel(consess)
            %        tmap_fn = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', k) '.nii']);
            tmap_clusters_fn = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', k) '_binary_clusters.nii']);
            %        [ptmap, ~, ~, ~] = fmrwhy_util_readOrientNifti(tmap_fn);
            nii = nii_tool('load', tmap_clusters_fn);
            func_roi_img = fmrwhy_util_createBinaryImg(nii.img, 0.1);

            overlap_img = anat_roi_img & func_roi_img;
            size_anat_roi = sum(anat_roi_img(:));
            size_func_roi = sum(func_roi_img(:));
            size_overlap = sum(overlap_img(:));
            dice_coeff = 2 * size_overlap / (size_anat_roi + size_func_roi);

            no_scaling = 1;
            str = consess{k}.tcon.name;
            saveAs_overlap_fn = fullfile(run_dir_stats, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-' str 'Overlaps' options.roi.(task).desc{j} '_roi.nii']);
            fmrwhy_util_saveNifti(saveAs_overlap_fn, double(overlap_img), tmap_clusters_fn, no_scaling);
            [poverlap, ~, ~, ~] = fmrwhy_util_readOrientNifti(saveAs_overlap_fn);
            plot_overlap_img = double(poverlap.nii.img);
            saveAs_fn = strrep(saveAs_overlap_fn, '.nii', '.png');
            roi_rgbcolors = [148, 239, 255];

            overlapmontage = fmrwhy_util_createStatsOverlayMontage(background_img, [], plot_overlap_img, 9, 1, '', 'gray', 'off', 'max', [], [], roi_rgbcolors, false, saveAs_fn);

        end
    end
