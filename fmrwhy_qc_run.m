function fmrwhy_qc_run(bids_dir, sub, ses, task, run, echo, options)

% Software dependences:
% - Matlab vX
% - SPM12 r7771
% - (TAPAS PhysIO Toolbox vX)

% Data required in order for this function to run as intended:
% - T1w (raw data)
% - Single run of fMRI timeseries (raw data)
% - Head motion/movement parameters derived from unprocessed data (i.e. from realignment of raw fMRI timeseries)
% - (T1w coregistered to template functional volume (from fmrwhy_preproc_structFunc.m))
% - (Segmentations of GM, WM, CSF in template functional volume space (from fmrwhy_preproc_structFunc.m))

% Code steps:
% 1. Get BIDS directory of run
% 2. Create relevant derivative directories (qc and preproc) if it doesn't exist yet
% 3. QC steps for fmri runs in subject folder:
%   - Head motion/movement parameters derived from unprocessed data (i.e. from realignment of raw fMRI timeseries)
%   - Framewise displacement
%   - Statistical measures / images (tsnr, variance, std, psc, DVARS)
%   - The plot
%   - PhysIO quality metrics and plots

% INPUT
% - flag to generate pictures or not
% - flag to generate nii images or not


% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);


% ------------------------
% SECTION A: ANATOMICAL QC
% ------------------------

% -------
% STEP 1: Contours of tissue masks on mean EPI / template EPI (/ anatomical image?)
% -------
% TODO: this should hide figures and only print them to png. currently figures are popping up.

mask_montage_fns = {'_GM_mask_montage', '_WM_mask_montage', '_CSF_mask_montage', '_brain_mask_montage'};
for i = 1:numel(mask_montage_fns)
    mask_montage_fns{i} = fullfile(options.anat_dir_qc, ['sub-' sub mask_montage_fns{i} '.png']);
end
run_montage = 0;
for i = 1:numel(mask_montage_fns)
    if ~exist(mask_montage_fns{i}, 'file')
        disp(['Mask overlay montage does not exist yet: ' mask_montage_fns{i}]);
        run_montage = 1;
    end
end
if run_montage
    disp('Creating mask overlay montages')
    montage_data = fmrwhy_qc_createMaskMontages(bids_dir, sub, 1);
    disp('Complete!')
    disp('---')
else
    disp('All mask overlay montages exist.')
    disp('---')
end

% -------
% STEP 2: Contours of anatomical ROIs on mean EPI / template EPI (/ anatomical image?)
% -------
[p1, frm1, rg1, dim1] = fmrwhy_util_readOrientNifti(options.rcoregest_anatomical_fn);
% Loop through all tasks in BIDS structure
for i = 1:numel(options.tasks)
    % Ignore the 'rest' task (assume there is no task ROI for this; have to change in future if RSnetworks available to be normalised or something)
    if strcmp(options.tasks{i}, 'rest') ~= 1
        % Loop through all ROIs for the particular task
        for j = 1:numel(options.roi.(options.tasks{i}).orig_fn)
            [p2, frm2, rg2, dim2] = fmrwhy_util_readOrientNifti(options.roi.(options.tasks{i}).rroi_fn{j});
            overlay_img = fmrwhy_util_createBinaryImg(p2.nii.img, 0.1);
            title = options.roi.(options.tasks{i}).name{j}
            saveAs_fn = fullfile(options.anat_dir_qc, ['sub-' sub '_space-individual_desc-' options.roi.(options.tasks{i}).desc{j} '_roi_montage.png']);
            if ~exist(saveAs_fn, 'file')
                fmrwhy_util_createOverlayMontage(p1.nii.img, overlay_img, 9, 1, title, 'gray', 'off', 'max', [], [255,0,0], saveAs_fn);
            else
                disp(['File already exists: ' saveAs_fn])
            end
        end
    end
end


% ------------------------
% SECTION B: FUNCTIONAL QC
% ------------------------

% -------
% STEP 1: Grab multiple regressors
% -------
%motion_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_motion.tsv']);
%framewise_displacement_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_fd.tsv']);
%tissue_regr_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_tissue.tsv']);
%physio_regr_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_physio.tsv']);
%confounds_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_regressors.tsv']);

% -------
% STEP 2: Calculate and generate statistical measures / images (tsnr, variance, std, psc, DVARS)
% -------
stats = fmrwhy_qc_createStatsOutput(bids_dir, sub, ses, task, run, echo, options);

% -------
% STEP 3: Create The Plot for whole brain
% -------
theplot_fn = {};
theplot_fn{1} = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-RO_grayplot.png']);
theplot_fn{2} = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-GSO_grayplot.png']);
if ~exist(theplot_fn{1}, 'file') || ~exist(theplot_fn{2}, 'file')
    fmrwhy_qc_createThePlot(bids_dir, sub, ses, task, run, echo, options);
end


% -------
% STEP 4: Create The Plot for rois
% -------

% Update workflow params with subject functional derivative filenames
options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);
% Ignore the 'rest' task (assume there is no task ROI for this; have to change in future if RSnetworks available to be normalised or something)
if strcmp(task, 'rest') ~= 1
    % Loop through all ROIs for the particular task
    for j = 1:numel(options.roi.(task).orig_fn)
        functional_fn = options.sfunctional_fn;
        desc = options.roi.(task).desc{j};
        saveAs_fn = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-' desc '_grayplot.png']);

        task_info.TR = options.firstlevel.(task).(['run' run]).sess_params.timing_RT;
        task_info.onsets = options.firstlevel.(task).(['run' run]).plot_params.cond_onset;
        task_info.durations = options.firstlevel.(task).(['run' run]).plot_params.cond_duration;
        task_info.precision = 1;

        if ~exist(saveAs_fn, 'file')
            trace_info = [];
            fmrwhy_util_thePlotROI(functional_fn, options.brain_mask_fn, options.roi.(task).rroi_fn{j}, task_info, trace_info, saveAs_fn)
        else
            disp(['File already exists: ' saveAs_fn])
        end
    end
else
    disp('---')
    disp('Not creating ROI timeseries plots for task = rest.')
    disp('---')
end