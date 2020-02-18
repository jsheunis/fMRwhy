% Pipeline to run through structural-functional preprocessing for a subjects

% Load/create required defaults
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
template_task = 'rest';
template_run = 1;
template_echo = 2;

% Input params
bids_dir = '/Volumes/Stephan_WD/NEUFEPME_data_BIDS';
sub = '001';
%ses = '';
%task = '';
%run = 1;
%echo = 2;
% BIDS structure values
BIDS = spm_BIDS(bids_dir);
subjects = spm_BIDS(BIDS,'subjects');
sessions = spm_BIDS(BIDS,'sessions');
runs = spm_BIDS(BIDS,'runs');
tasks = spm_BIDS(BIDS,'tasks');
types = spm_BIDS(BIDS,'types');
modalities = spm_BIDS(BIDS,'modalities');

%%
% Step 0.1: create necessary directory structure and copy required files
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

%%
% Step 0.2: create template functional volume nifti
template_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' num2str(template_run) '_space-individual_bold.nii'])

files = spm_BIDS(BIDS,'data', 'sub', sub ,'task', template_task, 'run', num2str(template_run), 'echo', num2str(template_echo),'type','bold')
if numel(files) ~= 1
    warning(['Expected a cell array with a single file, founc a cell array with ' numel(files) ' files.'])
else
    functional0_fn = [files{1} ',1'];
end
functional0_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' num2str(template_run) '_echo-' num2str(template_echo) '_bold.nii,1'])

fmrwhy_util_saveNifti(template_fn, spm_read_vols(spm_vol(functional0_fn)), functional0_fn, 'Template functional volume', 0)

%%
% Step 1: structural-functional-preproc:    - rtme_preproc_structFunc.m
fmrwhy_preproc_structFunc(bids_dir, sub, ses, template_task, template_run, template_echo);

%%
% Step 2: basic-functional-preproc:    - rtme_preproc_basicFunc.m
fmrwhy_preproc_basicFunc(sub, defaults);

% Step 3: anatomical-localizer-preproc:     - rtme_preproc_anatLocaliser.m
fmrwhy_preproc_anatLocaliser(sub, defaults);

% Step 3: functional-localizer-preproc:     - rtme_preproc_funcLocaliser.m
%                                           - rtme_preproc_generateRegressors.m
%                                           - rtme_preproc_generateRetroicor.m
%                                           - rtme_preproc_generateFDregr.m
%                                           - rtme_preproc_generateTissueSignals.m

% Step X: quality-preproc - rtme_preproc_qualityControl.m


for t = 1:numel(defaults.tasks)
%    disp(['Performing 3D volume realignment for: ' sub '_task-' tasks(t) '_run-' num2str(template_run)])
    rtme_preproc_funcLocaliser(sub, task, template_run, template_echo, defaults)
end


% Step 4: calculate-prior-measures-preproc - rtme_preproc_estimateParams.m
rtme_preproc_estimateParams(sub, defaults);










%### Pre-processing: peripheral data (RUN 1)
%
%1. Generate RETROICOR regressors from cardiac and respiratory traces of both runs (run 2 data to be used later) - PhysIO + Matlab
%
%
%### Pre-processing: functional (RUN 1)
%
%1. Task region localisation (using only middle echo [TE=28ms] timeseries):
%    1. Slice time correction
%    2. 3D volume realignment
%    3. Calculate framewise displacement from realignment params, select outliers using FD threshold (*which value or percentage?*)
%    4. Gaussian kernel smoothing (2*voxel size?)
%    5. GLM analysis incl:
%        1. AR(1) autoregressive filtering
%        2. Drift removal / high-pass (SPM cosine basis set)
%        3. Realignment params [+expansion?]
%        4. RETROICOR (+HRV, RTV?)
%        5. FD outlier binary regressor
%        6. *(global or tissue compartment signals???)*
%    6. Select t-stat peak within anatomically bound mask (from anatomy toolbox ROI)
%    7. Select N amount of voxels neighbouring peak voxel ==> ROI for real-time use
%
%2. T2*, S0, tSNR calculation from `run1_BOLD_rest` dataset (*is this sensible, as opposed to using RUN 1 task data?*):
%    1. Slice time correction on all three echo timeseries
%    2. 3D volume realignment on middle echo timeseries
%    3. Apply rigid body transformations from middle echo realignment parameters to echo 1 and echo 3 timeseries
%    6. T2* and S0 estimation (*check steps of tedana*):
%        1. *How to mask?*
%        2. Calculate timeseries average
%        3. Estimate T2* and S0 using log-linear fit of mono-exponential decay model
%        4. *Threshold?*
%    4. *Drift removal?*
%    5. tSNR calculation:
%        1. *How to mask?*
%        2. Mean / stddev












% Now, this step is executed once all settings are completed. We first
% check if the data has been preprocessed already. If so, we just load and
% name the variables. If not we run the standard preprocesing pipeline.
[d, fn, ext] = fileparts(structural_fn);
% We check if the data has been preprocessed by searching for the filename
% of one of the files that are generated during preprocessing (here, the
% resliced grey matter segmentation image)
if exist([d filesep 'rc1' fn ext], 'file')
    % Just load file/variable names, don't redo preprocessing
    disp('Preprocessing already done - loading variables')
    preproc_data = struct;
    [d, fn, ext] = fileparts(structural_fn);
    preproc_data.forward_transformation = [d filesep 'y_' fn ext];
    preproc_data.inverse_transformation = [d filesep 'iy_' fn ext];
    preproc_data.gm_fn = [d filesep 'c1' fn ext];
    preproc_data.wm_fn = [d filesep 'c2' fn ext];
    preproc_data.csf_fn = [d filesep 'c3' fn ext];
    preproc_data.bone_fn = [d filesep 'c4' fn ext];
    preproc_data.soft_fn = [d filesep 'c5' fn ext];
    preproc_data.air_fn = [d filesep 'c6' fn ext];
    preproc_data.rstructural_fn = [d filesep 'r' fn ext];
    preproc_data.rgm_fn = [d filesep 'rc1' fn ext];
    preproc_data.rwm_fn = [d filesep 'rc2' fn ext];
    preproc_data.rcsf_fn = [d filesep 'rc3' fn ext];
    preproc_data.rbone_fn = [d filesep 'rc4' fn ext];
    preproc_data.rsoft_fn = [d filesep 'rc5' fn ext];
    preproc_data.rair_fn = [d filesep 'rc6' fn ext];

    % Check if ROIs are specified to be in native space and run warping or
    % not based on setting
    if ROI_native
        % ROIs are already in native space, no warping necessary
        preproc_data.ROI_fns = ROI_fns;
    else
        % This part was hardcoded for testing purposes, it doesnt actually
        % call the warping functionality here, which it should do. TODO
        for roi = 1:(N_ROIs-N_RNOIs)
            [droi, fnroi, extroi] = fileparts(ROI_fns{roi});
            preproc_data.wROI_fns{roi} = [droi filesep 'w' fnroi extroi];
            preproc_data.rwROI_fns{roi} = [droi filesep 'rw' fnroi extroi];
        end
        preproc_data.ROI_fns = preproc_data.rwROI_fns;
    end
else
    % If preproc not done, call preprocessing script
    preproc_data = onlineBrain_preRtPreProc(functional0_fn, structural_fn, spm_dir);

    % Also do ROI warping to get everything in same space
    if ROI_native
        % ROIs are already in native space, no warping necessary
        preproc_data.ROI_fns = ROI_fns;
    else
        % ROIs are in MNI space, warping and reslicing necessary
        % Warp MNI space rois to functional space,...
        spm_normalizeWrite_jsh(preproc_data.inverse_transformation, ROI_fns(1:(end-N_RNOIs)));
        for roi = 1:(N_ROIs-N_RNOIs)
            [droi, fnroi, extroi] = fileparts(ROI_fns{roi});
            preproc_data.wROI_fns{roi} = [droi filesep 'w' fnroi extroi];
        end
        % ... then reslice
        spm('defaults','fmri');
        spm_jobman('initcfg');
        reslice = struct;
        % Ref
        reslice.matlabbatch{1}.spm.spatial.coreg.write.ref = {functional0_fn};
        % Source
        source_fns = preproc_data.wROI_fns;
        reslice.matlabbatch{1}.spm.spatial.coreg.write.source = source_fns';
        % Roptions
        reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
        reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
        reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
        reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
        % Run
        spm_jobman('run',reslice.matlabbatch);
        for roi = 1:(N_ROIs-N_RNOIs)
            [droi, fnroi, extroi] = fileparts(ROI_fns{roi});
            preproc_data.rwROI_fns{roi} = [droi filesep 'rw' fnroi extroi];
        end
        preproc_data.ROI_fns = preproc_data.rwROI_fns;
    end

end