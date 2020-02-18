% A custom workflow that does structFunc and basicFunc preprocessing and QC fof a single subject in the NEUFEP study

% Code steps:
% 1. Define template/default variables, directories and filenames
% 2. Specify


%--------------------------------------------------------------------------

% Load/create required defaults
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
bids_dir = '';
sub = '001';
ses = '';
template_task = 'rest';
template_run = '1';
template_echo = '2';

% BIDS structure values
BIDS = spm_BIDS(bids_dir);
subjects = spm_BIDS(BIDS,'subjects');
sessions = spm_BIDS(BIDS,'sessions');
runs = spm_BIDS(BIDS,'runs');
tasks = spm_BIDS(BIDS,'tasks');
types = spm_BIDS(BIDS,'types');
modalities = spm_BIDS(BIDS,'modalities');

% Directory and content setup
deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
qc_dir = fullfile(deriv_dir, 'fmrwhy-qc');
sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);
sub_dir_qc = fullfile(qc_dir, ['sub-' sub]);
if ~exist(sub_dir_preproc, 'dir')
    mkdir(sub_dir_preproc)
    sub_dir_BIDS = fullfile(bids_dir, ['sub-' sub]);
    copyfile(sub_dir_BIDS, sub_dir_preproc)
end
if ~exist(sub_dir_qc, 'dir')
    mkdir(sub_dir)
end

% Get anatomical file
anatomical_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_T1w.nii']);

% -------
% Get structfunc preproc filenames
% -------
% Outputs after coregistering T1w image
coregest_anatomical_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-coregEst_T1w.nii']);
% Outputs after segmenting coregistered T1w image
gm_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-GM_probseg.nii']);
wm_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-WM_probseg.nii']);
csf_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-CSF_probseg.nii']);
bone_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-bone_probseg.nii']);
soft_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-soft_probseg.nii']);
air_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-air_probseg.nii']);
indiv_to_mni_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_desc-IndivToMNI_transform.nii']); % forward transform
mni_to_indiv_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_desc-MNItoIndiv_transform.nii']); % inverse transform
probseg_fns = {gm_probseg_fn, wm_probseg_fn, csf_probseg_fn, bone_probseg_fn, soft_probseg_fn, air_probseg_fn};
transform_fns = {indiv_to_mni_fn, mni_to_indiv_fn};
% Outputs after reslicing segments and coregistered T1w image
rcoregest_anatomical_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-coregEstResl_T1w.nii']);
rgm_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rGM_probseg.nii']);
rwm_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rWM_probseg.nii']);
rcsf_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rCSF_probseg.nii']);
rbone_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rbone_probseg.nii']);
rsoft_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rsoft_probseg.nii']);
rair_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rair_probseg.nii']);
rprobseg_fns = {rgm_probseg_fn, rwm_probseg_fn, rcsf_probseg_fn, rbone_probseg_fn, rsoft_probseg_fn, rair_probseg_fn};
rall_fns = {rcoregest_anatomical_fn, rgm_probseg_fn, rwm_probseg_fn, rcsf_probseg_fn, rbone_probseg_fn, rsoft_probseg_fn, rair_probseg_fn};
% Outputs after creating masks
gm_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-GM_mask.nii']);
wm_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-WM_mask.nii']);
csf_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-CSF_mask.nii']);
brain_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-brain_mask.nii']);
mask_fns = {gm_mask_fn, wm_mask_fn, csf_mask_fn, brain_mask_fn};
% -------
% Get basicfunc preproc filenames
% -------
motion_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_desc-confounds_motion.tsv']);
afunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-apreproc_bold.nii']);
rfunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rpreproc_bold.nii']);
rafunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_bold.nii']);
sfunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-spreproc_bold.nii']);
srfunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srpreproc_bold.nii']);
srafunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srapreproc_bold.nii']);
framewise_displacement_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_fd.tsv']);
tissue_regr_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_tissue.tsv']);
physio_regr_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_physio.tsv']);
confounds_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_regressors.tsv']);


% -------
% STEP 0 -- Create functional template
% -------
% Define template filename
template_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_space-individual_bold.nii'])
% Create, if it does not exist
if ~exist(template_fn, 'file')
    functional0_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_echo-' template_echo '_bold.nii,1'])
    fmrwhy_util_saveNifti(template_fn, spm_read_vols(spm_vol(functional0_fn)), functional0_fn, 'Template functional volume', 0)
end

% -------
% STEP 1 -- Structural-functional preprocessing: fmrwhy_preproc_structFunc.m
% -------
struct_func_out_fns = {coregest_anatomical_fn, probseg_fns, transform_fns, rall_fns, mask_fns};
run_structFunc = 0;
for i = 1:numel(struct_func_out_fns)
    if ~exist(struct_func_out_fns{i}, 'file')
        disp(['Structural-funcional preprocessing output file does not exist yet: ' struct_func_out_fns{i}]);
        run_structFunc = 1;
    end
end
if run_structFunc
    fmrwhy_preproc_structFunc(bids_dir, sub, ses, template_task, template_run, template_echo);
    disp('Complete!')
    disp('---')
else
    disp('Structural-funcional preprocessing already completed.')
    disp('---')
end

% -------
% STEP 2 -- Basic functional preprocessing: fmrwhy_preproc_basicFunc.m
% -------
basic_func_out_fns = {motion_fn, afunctional_fn, rfunctional_fn, rafunctional_fn, sfunctional_fn, srfunctional_fn, srafunctional_fn};
run_basicFunc = 0;
for i = 1:numel(basic_func_out_fns)
    if ~exist(basic_func_out_fns{i}, 'file')
        disp(['Basic funcional preprocessing output file does not exist yet: ' basic_func_out_fns{i}]);
        run_basicFunc = 1;
    end
end
if run_basicFunc
    fmrwhy_preproc_basicFunc(bids_dir, sub, ses, template_task, template_run, template_echo);
    disp('Complete!')
    disp('---')
else
    disp('Basic funcional preprocessing already completed.')
    disp('---')
end


% -------
% STEP 3 -- Quality control pipeline: fmrwhy_qc_.m
% -------
% Step 2: quality-preproc - fmrwhy_qc_sub.m
qc_out_fns; % some file that is generated by the qc procedure (subject level or run level?)
run_qc = 0;
for i = 1:numel(qc_out_fns)
    if ~exist(qc_out_fns{i}, 'file')
        disp(['Basic funcional preprocessing output file does not exist yet: ' qc_out_fns{i}]);
        run_qc = 1;
    end
end
if run_qc
    fmrwhy_qc_run(bids_dir, sub, ses, template_task, template_run, template_echo);
    disp('Complete!')
    disp('---')
else
    disp('Basic funcional preprocessing already completed.')
    disp('---')
end


%
%% Step 3: anatomical-localizer-preproc:     - rtme_preproc_anatLocaliser.m
%fmrwhy_preproc_anatLocaliser(sub, defaults);

% Step 3: functional-localizer-preproc:     - rtme_preproc_funcLocaliser.m
%                                           - rtme_preproc_generateRegressors.m
%                                           - rtme_preproc_generateRetroicor.m
%                                           - rtme_preproc_generateFDregr.m
%                                           - rtme_preproc_generateTissueSignals.m



%
%for t = 1:numel(defaults.tasks)
%%    disp(['Performing 3D volume realignment for: ' sub '_task-' tasks(t) '_run-' template_run])
%    rtme_preproc_funcLocaliser(sub, task, template_run, template_echo, defaults)
%end


% Step 4: calculate-prior-measures-preproc - rtme_preproc_estimateParams.m
%rtme_preproc_estimateParams(sub, defaults);










%%### Pre-processing: peripheral data (RUN 1)
%%
%%1. Generate RETROICOR regressors from cardiac and respiratory traces of both runs (run 2 data to be used later) - PhysIO + Matlab
%%
%%
%%### Pre-processing: functional (RUN 1)
%%
%%1. Task region localisation (using only middle echo [TE=28ms] timeseries):
%%    1. Slice time correction
%%    2. 3D volume realignment
%%    3. Calculate framewise displacement from realignment params, select outliers using FD threshold (*which value or percentage?*)
%%    4. Gaussian kernel smoothing (2*voxel size?)
%%    5. GLM analysis incl:
%%        1. AR(1) autoregressive filtering
%%        2. Drift removal / high-pass (SPM cosine basis set)
%%        3. Realignment params [+expansion?]
%%        4. RETROICOR (+HRV, RTV?)
%%        5. FD outlier binary regressor
%%        6. *(global or tissue compartment signals???)*
%%    6. Select t-stat peak within anatomically bound mask (from anatomy toolbox ROI)
%%    7. Select N amount of voxels neighbouring peak voxel ==> ROI for real-time use
%%
%%2. T2*, S0, tSNR calculation from `run1_BOLD_rest` dataset (*is this sensible, as opposed to using RUN 1 task data?*):
%%    1. Slice time correction on all three echo timeseries
%%    2. 3D volume realignment on middle echo timeseries
%%    3. Apply rigid body transformations from middle echo realignment parameters to echo 1 and echo 3 timeseries
%%    6. T2* and S0 estimation (*check steps of tedana*):
%%        1. *How to mask?*
%%        2. Calculate timeseries average
%%        3. Estimate T2* and S0 using log-linear fit of mono-exponential decay model
%%        4. *Threshold?*
%%    4. *Drift removal?*
%%    5. tSNR calculation:
%%        1. *How to mask?*
%%        2. Mean / stddev
%
%
%
%
%
%
%
%
%
%
%
%
%% Now, this step is executed once all settings are completed. We first
%% check if the data has been preprocessed already. If so, we just load and
%% name the variables. If not we run the standard preprocesing pipeline.
%[d, fn, ext] = fileparts(structural_fn);
%% We check if the data has been preprocessed by searching for the filename
%% of one of the files that are generated during preprocessing (here, the
%% resliced grey matter segmentation image)
%if exist([d filesep 'rc1' fn ext], 'file')
%    % Just load file/variable names, don't redo preprocessing
%    disp('Preprocessing already done - loading variables')
%    preproc_data = struct;
%    [d, fn, ext] = fileparts(structural_fn);
%    preproc_data.forward_transformation = [d filesep 'y_' fn ext];
%    preproc_data.inverse_transformation = [d filesep 'iy_' fn ext];
%    preproc_data.gm_fn = [d filesep 'c1' fn ext];
%    preproc_data.wm_fn = [d filesep 'c2' fn ext];
%    preproc_data.csf_fn = [d filesep 'c3' fn ext];
%    preproc_data.bone_fn = [d filesep 'c4' fn ext];
%    preproc_data.soft_fn = [d filesep 'c5' fn ext];
%    preproc_data.air_fn = [d filesep 'c6' fn ext];
%    preproc_data.rstructural_fn = [d filesep 'r' fn ext];
%    preproc_data.rgm_fn = [d filesep 'rc1' fn ext];
%    preproc_data.rwm_fn = [d filesep 'rc2' fn ext];
%    preproc_data.rcsf_fn = [d filesep 'rc3' fn ext];
%    preproc_data.rbone_fn = [d filesep 'rc4' fn ext];
%    preproc_data.rsoft_fn = [d filesep 'rc5' fn ext];
%    preproc_data.rair_fn = [d filesep 'rc6' fn ext];
%
%    % Check if ROIs are specified to be in native space and run warping or
%    % not based on setting
%    if ROI_native
%        % ROIs are already in native space, no warping necessary
%        preproc_data.ROI_fns = ROI_fns;
%    else
%        % This part was hardcoded for testing purposes, it doesnt actually
%        % call the warping functionality here, which it should do. TODO
%        for roi = 1:(N_ROIs-N_RNOIs)
%            [droi, fnroi, extroi] = fileparts(ROI_fns{roi});
%            preproc_data.wROI_fns{roi} = [droi filesep 'w' fnroi extroi];
%            preproc_data.rwROI_fns{roi} = [droi filesep 'rw' fnroi extroi];
%        end
%        preproc_data.ROI_fns = preproc_data.rwROI_fns;
%    end
%else
%    % If preproc not done, call preprocessing script
%    preproc_data = onlineBrain_preRtPreProc(functional0_fn, structural_fn, spm_dir);
%
%    % Also do ROI warping to get everything in same space
%    if ROI_native
%        % ROIs are already in native space, no warping necessary
%        preproc_data.ROI_fns = ROI_fns;
%    else
%        % ROIs are in MNI space, warping and reslicing necessary
%        % Warp MNI space rois to functional space,...
%        spm_normalizeWrite_jsh(preproc_data.inverse_transformation, ROI_fns(1:(end-N_RNOIs)));
%        for roi = 1:(N_ROIs-N_RNOIs)
%            [droi, fnroi, extroi] = fileparts(ROI_fns{roi});
%            preproc_data.wROI_fns{roi} = [droi filesep 'w' fnroi extroi];
%        end
%        % ... then reslice
%        spm('defaults','fmri');
%        spm_jobman('initcfg');
%        reslice = struct;
%        % Ref
%        reslice.matlabbatch{1}.spm.spatial.coreg.write.ref = {functional0_fn};
%        % Source
%        source_fns = preproc_data.wROI_fns;
%        reslice.matlabbatch{1}.spm.spatial.coreg.write.source = source_fns';
%        % Roptions
%        reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
%        reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
%        reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
%        reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
%        % Run
%        spm_jobman('run',reslice.matlabbatch);
%        for roi = 1:(N_ROIs-N_RNOIs)
%            [droi, fnroi, extroi] = fileparts(ROI_fns{roi});
%            preproc_data.rwROI_fns{roi} = [droi filesep 'rw' fnroi extroi];
%        end
%        preproc_data.ROI_fns = preproc_data.rwROI_fns;
%    end
%
%end