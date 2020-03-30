%%
% This script initialises the settings necessary for:
% ---------------------- %
% File/directory locations
% Data initialisation
% Processing and execution settings
% Experimental and Protocol parameters
% Data setup (scripts)
% Pre-real-time preprocessing (scripts)
% ---------------------- %
% It also runs some code/scripts to calculate variables necessary for
% online use from the initialised data/variables.


%% File/directory locations

matlab_dir = '/Users/jheunis/Documents/MATLAB'; % Specify MATLAB code directory
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12'; % Specify SPM installation directory

%% Data setup

% % --- For NEUFEP-ME data (single echo test) --- %
data_dir = '/Users/jheunis/Desktop/sample-data/sub-neufepmetest'; % Specify parent directory that contains all data
sub = 'sub-pilot';
sub_dir = [data_dir filesep sub]; % Specify specific subject directory
func_dir = [sub_dir filesep 'func'];
anat_dir = [sub_dir filesep 'anat']; 
roi_dir = [sub_dir filesep 'roi'];
functional0_fn = [func_dir filesep sub '_task-rest_run-1_echo-2.nii,1']; % Functional template from pre-real-time scan (middle scan from multi-echo*3)
structural_fn = [anat_dir filesep sub '_T1w.nii']; % Structural scan, from pre-real-time

ROI_native = true; % are ROIs already in native space with matching resolution?
N_RNOIs = 4; % Reference ROIs: 1-GM, 2-WM, 3-CSF, 4-masked brain (=GM+WM+CSF), 5-background slice/region
N_ROIs = 2+N_RNOIs; % number of ROIs supplied + N_RNOIs
ROI_fns = cell(1,N_ROIs);
% ROI_fns{1} = [roi_dir filesep 'rwBilateral_Amygdala_allregions.nii'];
% ROI_fns{1} = [roi_dir filesep 'rwLeft_Amygdala_allregions.nii'];
% ROI_fns{2} = [roi_dir filesep 'rwRight_Amygdala_allregions.nii'];
ROI_fns{1} = [roi_dir filesep 'rwLeft_Motor_4a_4p.nii'];
ROI_fns{2} = [roi_dir filesep 'rwRight_Motor_4a_4p.nii'];
ROI_names = cell(1,N_ROIs);
% ROI_names{1} = 'Bilateral amygdala';
% ROI_names{1} = 'Left amygdala';
% ROI_names{2} = 'Right amygdala';
% ROI_names{1} = 'Left motor';
% ROI_names{2} = 'Right motor';

Ne = 1; % number of EPI echoes

% Eref = 2; % reference echo, for realignment
% TE = [12 35 58]; % Echo times in ms
% use_echo = 2; % specify if only a single echo should be used for real-time
% template_echo = 2; % echo to be used for motion parameters, for ME
% T2star_thresh = 100; % threshold for maximum T2star after estimation
voxel_size = [3.5 3.5 3.5];
smoothing_kernel    = [7 7 7];
Nt =   210; % NrOfVolumes % VolumesNumber
N_skip = 0; % nrSkipVol
N_start = N_skip + 1;
Ndyn = Nt - N_skip;
TR = 2;
TR_skip = 2; % amount of TRs to skip at start of baseline block to exclude effect of HRF on baseline mean calculation
NF_cond = 2; % location of task/nf condition in design matrix (SPM structure)
lCond = 2; % number of conditions (should actually be derived form SPM mat)
% Emotion run 1
% task_onsets = [13; 33; 53; 73; 93; 113; 133; 153; 173; 193]; % onsets (in scan number) of task blocks in experimental paradigm
% task_durations = [8; 8; 8; 8; 8; 8; 8; 8; 8; 8];
% baseline_onsets = [1; 21; 41; 61; 81; 101; 121; 141; 161; 181; 201];
% baseline_durations = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10];
% Emotion run 2, and motor
task_onsets = [11; 31; 51; 71; 91; 111; 131; 151; 171; 191]; % onsets (in scan number) of task blocks in experimental paradigm
task_durations = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10];
baseline_onsets = [1; 21; 41; 61; 81; 101; 121; 141; 161; 181; 201];
baseline_durations = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10];

Nslice = 17;
rotateDir = [0 0 1];
rotateVal = 1;
rotateDeg = 270;
showMontage = true;

% % % --- For ME test 23 July --- %
% data_dir = 'C:\Users\HeunisS\Desktop\Matlabdata'; % Specify parent directory that contains all data
% dump_dir = 'C:\Drindumps'; % Specify parent directory where real-time nifti files are "dumped"
% sub = 'test_25Jul2019';
% sub_dir = [data_dir filesep sub]; % Specify specific subject directory
% if ~exist(sub_dir, 'dir')
%    mkdir(sub_dir)
% end
% func_dir = [sub_dir filesep 'func']; % Directory into which real-time images are copied
% if ~exist(func_dir, 'dir')
%    mkdir(func_dir)
% end
% functional0_fn = [sub_dir filesep '.nii']; % Functional template from pre-real-time scan (middle scan from multi-echo*3)
% structural_fn = [sub_dir filesep '.nii']; % Structural scan, from pre-real-time
% 
% ROI_native = false; % are ROIs already in native space with matching resolution?
% N_ROIs = 0+4; % number of ROIs supplied; last 4 are 1-GM, 2-WM, 3-CSF, 4-masked brain (=GM+WM+CSF; or other if chosen differently)
% ROI_fns = cell(1,N_ROIs);
% Ne = 3; % number of echoes per volume
% Eref = 2; % reference echo, for realignment
% TE = [12 35 58]; % Echo times in ms
% use_echo = 2; % specify if only a single echo should be used for real-time
% template_echo = 2; % echo to be used for motion parameters, for ME
% T2star_thresh = 100; % threshold for maximum T2star after estimation
% voxel_size = [3.5 3.5 4.5];
% smoothing_kernel    = [7 7 7];
% Nt =   208; % NrOfVolumes % VolumesNumber
% N_skip = 0; % nrSkipVol
% N_start = N_skip + 1;
% Ndyn = Nt - N_skip;
% TR = 2;
% TR_skip = 2; % amount of TRs to skip at start of baseline block to exclude effect of HRF on baseline mean calculation
% NF_cond = 2; % location of task/nf condition in design matrix (SPM structure)
% lCond = 2; % number of conditions (should actually be derived form SPM mat)
% task_onsets = [17; 49; 81; 113; 145; 177];
% task_durations = [16; 16; 16; 16; 16; 16];
% baseline_onsets = [1; 33; 65; 97; 129; 161; 193];
% baseline_durations = [16; 16; 16; 16; 16; 16; 16];
% Nslice = 17;
% rotateDir = [0 0 1];
% rotateVal = 1;
% rotateDeg = 270;

% % --- For OpenNFT sample data --- %
% data_dir = '/Users/jheunis/Desktop/All/Code and data tests/neu3carttest'; % Specify parent directory that contains all data
% sub_dir = [data_dir filesep 'sub-opennft']; % Specify specific subject directory
% functional0_fn      =   [sub_dir filesep 'template_func.nii']; % Functional template from pre-real-time scan
% structural_fn = [sub_dir filesep 'structScan_PSC.nii']; % Structural scan, from pre-real-time
% ROI_native = true; % are ROIs already in native space with matching resolution?
% N_RNOIs = 5; % Reference ROIs: 1-GM, 2-WM, 3-CSF, 4-masked brain (=GM+WM+CSF), 5-background slice/region
% N_ROIs = 2+N_RNOIs; % number of ROIs supplied + N_RNOIs
% ROI_fns = cell(1,N_ROIs);
% ROI_fns{1} = [sub_dir filesep 'lROI_1.nii'];
% ROI_fns{2} = [sub_dir filesep 'rROI_2.nii'];
% Ne = 1; % number of EPI echoes
% voxel_size = [2.973 2.973 3.75];
% smoothing_kernel    = [6 6 6];
% Nt =   155; % NrOfVolumes % VolumesNumber
% N_skip = 5; % nrSkipVol
% N_start = N_skip + 1;
% Ndyn = Nt - N_skip;
% TR = 2;
% TR_skip = 2; % amount of TRs to skip at start of baseline block to exclude effect of HRF on baseline mean calculation
% NF_cond = 2; % location of task/nf condition in design matrix (in SPM.mat structure)
% lCond = 2; % number of conditions (should actually be derived form SPM.mat)
% timing_units = 'scans';
% task_onsets = [11; 31; 51; 71; 91; 111; 131]; % onsets (in scan number) of task blocks in experimental paradigm
% task_durations = [10; 10; 10; 10; 10; 10; 10];
% baseline_onsets = [1; 21; 41; 61; 81; 101; 121; 141];
% baseline_durations = [10; 10; 10; 10; 10; 10; 10; 10];
% Nslice = 17;
% rotateDir = [0 0 1];
% rotateVal = 1;
% rotateDeg = 270;
% showMontage = true;


% % --- For LI Finger Tapping --- %
% data_dir            =   '/Users/jheunis/Desktop/All/Code and data tests/rt-fmri-dev-tests/LI_finger_tapping';
% sub = '02';
% sub_dir = [data_dir filesep 'sub-' sub];
% functional0_fn      =   [sub_dir filesep 'sub-' sub '_bold.nii,1']; % Functional template from pre-real-time scan
% structural_fn = [sub_dir filesep 'sub-' sub '_T1w.nii']; % Structural scan, from pre-real-time
% ROI_native = false; % are ROIs already in native space with matching resolution?
% N_ROIs = 4+4; % number of ROIs supplied; last 4 are 1-GM, 2-WM, 3-CSF, 4-masked brain (=GM+WM+CSF; or other if chosen differently)
% ROI_fns = cell(1,N_ROIs);
% ROI_fns{1} = [data_dir filesep 'Left_Motor_4a.nii'];
% ROI_fns{2} = [data_dir filesep 'Right_Motor_4a.nii'];
% ROI_fns{3} = [data_dir filesep 'Left_Motor_4p.nii'];
% ROI_fns{4} = [data_dir filesep 'Right_Motor_4p.nii'];
% Ne = 1; % number of EPI echoes
% voxel_size = [1.75 1.75 3];
% smoothing_kernel    = [5 5 5];
% Nt =   160; % NrOfVolumes % VolumesNumber
% N_skip = 0; % nrSkipVol
% N_start = N_skip + 1;
% Ndyn = Nt - N_skip;
% TR = 3;
% TR_skip = 2; % amount of TRs to skip at start of baseline block to exclude effect of HRF on baseline mean calculation
% NF_cond = 2; % location of task/nf condition in design matrix (SPM structure)
% lCond = 2; % number of conditions (should actually be derived form SPM mat)
% task_onsets = [11; 31; 51; 71; 91; 111; 131; 151];
% task_durations = [10; 10; 10; 10; 10; 10; 10; 10];
% baseline_onsets = [1; 21; 41; 61; 81; 101; 121; 141];
% baseline_durations = [10; 10; 10; 10; 10; 10; 10; 10];
% Nslice = 34;
% rotateDir = [0 0 1];
% rotateVal = -1;
% rotateDeg = 90;


%% Processing and execution settings
% (Note: many of these settings relate to processing steps that were
% implemented for this application based on OpenNFT code, so one might see
% the similarities when inspecting OpenNFT. For the new application in
% Python, similar settings might not be necessary, so it is not necessarily
% true that all lines of code in Matlab will form part of the eventual
% Python application.

% Processing steps to include: AR(1), cGLM, iGLM
isIGLM = false;
iglmAR1 = false;
isRegrIGLM = true;
isMotionRegr = true;
isHighPass = true;
isLinRegr = true;
isPhysRegr = false;
cglmAR1 = true;

% ?
nrBlocksInSlidingWindow = 100; % i.e disabled

% SPM motion correction settings
flagsSpmRealign = struct('quality',.9,'fwhm',5,'sep',4,...
    'interp',4,'wrap',[0 0 0],'rtm',0,'PW','','lkp',1:6);
flagsSpmReslice = struct('quality',.9,'fwhm',5,'sep',4,...
    'interp',4,'wrap',[0 0 0],'mask',1,'mean',0,'which', 2);

% SPM AR(1) settings
aAR1 = 0.2; % default SPM value

% Kalman filter settings
S = struct;
S.Q = 0;
S.P = S.Q;
S.x = 0;
fPositDerivSpike = 0;
fNegatDerivSpike = 0;
S(1:N_ROIs) = S;
fPositDerivSpike(1:N_ROIs) = fPositDerivSpike;
fNegatDerivSpike(1:N_ROIs) = fNegatDerivSpike;

% Scaling settings
fLockedTempl = 0; % 0 = update, 1 = fixed

tmp_posMin(1:N_ROIs) = 0;
tmp_posMax(1:N_ROIs) = 0;

% Empty variables init
rawTimeSeries = [];
kalmanProcTimeSeries = [];
displRawTimeSeries = [];
scalProcTimeSeries = [];
emaProcTimeSeries = [];
posMin = [];
posMax = [];
mposMax = [];
mposMin = [];
blockNF = 0;
firstNF = 0;

%% Protocol parameters (derived/calculated)

basBlockLength = task_durations(1); % Here hardcoded to be same as task block
% Logic for creating binary task-baseline vector, for determining if
% online iteration is part of task or baseline in terms of stimulus
baseline_mean = zeros(2,length(task_onsets)+1);
task_design = zeros(1, Ndyn);
for n = 1:length(task_onsets)
    for m = task_onsets(n):(task_onsets(n)+task_durations(n)-1)
        task_design(m) = 1;
    end
end
% Condition-encoded task-baseline vector
vectEncCond = task_design + 1; % 1 = index in baseline block; 2 = index in task block
% Get indices for task and baseline blocks
task_blocks = bwlabel(task_design);
baseline_blocks = bwlabel(~task_design);
baseline_design = (~task_design);
I_baseline = [];
for n = 1:length(baseline_onsets)
    I_block = find(baseline_blocks == n);
    if n > 1
        I_block = I_block((1+TR_skip):end);
    end
    I_baseline = [I_baseline I_block];
end
% Calculate HRF-convolved tas and basline designs. The HRF (hemodynamic
% response function) tells us that blood (and hence oxygen level) takes
% some seconds to respond to stimulus-driven neuronal activity. Thus when
% we e.g. see a stimulus, the oxygen level change due to the neuronal
% activity will only be measurable with fMRI a few seconds later. So we do
% convolution of the task design to have a more accurate expectatioon for
% the fMRI signal change.
hrf = spm_hrf(TR);
convolved_task_design = spm_Volterra(struct('u', task_design', 'name', {{'task'}}), hrf);
convolved_baseline_design = spm_Volterra(struct('u', double((~task_design)'), 'name', {{'baseline'}}), hrf);
tmpSpmDesign = convolved_task_design;
% Calculate cosine basis set for drift regressors. In fMRI analysis, GLM is
% used to model the data and we also include noise regressors in the model
% because we know fMRI data is very noisy. SPM contains some functionality
% that is used below to model noise regressors with varying frequancy
% content
K.HParam = 128;
K.RT = TR;
k    = Ndyn;
n    = fix(2*(k*K.RT)/K.HParam + 1);
X0   = spm_dctmtx(k,n);
K.X0 = X0(:,2:end);
cosine_basis_set = K.X0; % cosine basis set for drift regressors
% K = SPM.xX.K; % OPenNFT' way to access cosine basis set for drift regressors

% If AR1 selected for iGLM, filter the task and baseline regressors so that
% they don't have to be filtered in real-time
basFct = [convolved_baseline_design convolved_task_design];
if iglmAR1
    basFct = arRegr(aAR1, basFct);
end
nrBasFct = size(basFct,2);

% Other OpenNFT values without explanation, but necessary
spmMaskTh = 80*ones(Ndyn,1)'; %80*ones(size(SPM.xM.TH));% mean(SPM.xM.TH)*ones(size(SPM.xM.TH)); % SPM.xM.TH;
pVal = .01;
statMap3D_iGLM = [];
tContr = [0; 1]; % Baseline; Task]  TODO

%% Data setup

% Volume dimensions, and reference image
funcref_spm = spm_vol(functional0_fn);
funcref_3D  = spm_read_vols(funcref_spm);
[Nx, Ny, Nz] = size(funcref_3D);
Nvox = Nx*Ny*Nz;

% Get motion realignment template data and volume
dicomInfoVox   = sqrt(sum((funcref_spm.mat(1:3,1:3)).^2));
fwhm = smoothing_kernel ./ dicomInfoVox;
A0=[];x1=[];x2=[];x3=[];wt=[];deg=[];b=[];
% MULTI-ECHO TODO
R(1,1).mat = funcref_spm.mat;
R(1,1).dim = funcref_spm.dim;
R(1,1).Vol = funcref_3D;

% Regressor setup for iGLM/cGLM
% Constant and Linear regressors
constRegrFull = ones(Ndyn,1);
linRegrFull = zscore((1:double(Ndyn))');

% Initial data setup for iGLM
nrBasFctRegr = 1; % Constant regressor is always included
if isMotionRegr
    nrBasFctRegr = nrBasFctRegr + 6;
end
if isLinRegr
    nrBasFctRegr = nrBasFctRegr + 1;
end
if isHighPass
    nrBasFctRegr = nrBasFctRegr + size(cosine_basis_set,2);
end
Cn = zeros(nrBasFct + nrBasFctRegr);
Dn = zeros(Nvox, nrBasFct + nrBasFctRegr);
s2n = zeros(Nvox, 1);
tn = zeros(Nvox, 1);
tTh = zeros(Ndyn, 1);
dyntTh = 0;
statMapVect = zeros(Nvox, 1);
statMap3D = zeros(Nx, Ny, Nz);
statMap4D = cell(1,Ndyn);
tContr = [tContr; zeros(nrBasFctRegr,1)];

% AR(1) for cGLM in signal preproessing
if ~iglmAR1
    spmDesign = tmpSpmDesign;
else
    spmDesign = arRegr(aAR1, tmpSpmDesign);
end
mf = [];
npv = 0;
statMapCreated = 0;
nrRegrDesign = size(spmDesign,2);
% number of regressors of no interest to correct with cGLM
nrRegrToCorrect = 8; % 6 MC regressors, linear trend, constant


%% Pre-real-time preproc
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

% Whole brain masking from the segmentations that resulted from preproc
[GM_img_bin, WM_img_bin, CSF_img_bin] = onlineBrain_getSegments(preproc_data.rgm_fn, preproc_data.rwm_fn, preproc_data.rcsf_fn, 0.5);
I_GM = find(GM_img_bin);
I_WM = find(WM_img_bin);
I_CSF = find(CSF_img_bin);
mask_3D = GM_img_bin | WM_img_bin | CSF_img_bin;
I_full_mask = find(mask_3D);
I_mask = I_full_mask;
N_maskvox = numel(I_mask);

% ROI masking (all ROIs prespecified)
ROI_img = cell(1,N_ROIs);
for r = 1:(N_ROIs-N_RNOIs)
    ROI_img{r} = spm_read_vols(spm_vol(preproc_data.ROI_fns{r}));
    ROI_img{r}(isnan(ROI_img{r})) = 0;
    ROI_img{r} = (ROI_img{r} >= 0.1);
    I_roi{r} = find(ROI_img{r}(:) & GM_img_bin(:));
end
ROI_img{r+1} = GM_img_bin;
I_roi{r+1} = I_GM;
ROI_names{r+1} = 'GM';
ROI_img{r+2} = WM_img_bin;
I_roi{r+2} = I_WM;
ROI_names{r+2} = 'WM';
ROI_img{r+3} = CSF_img_bin;
I_roi{r+3} = I_CSF;
ROI_names{r+3} = 'CSF';
ROI_img{r+4} = mask_3D;
I_roi{r+4} = I_mask;
ROI_names{r+4} = 'Brain';
% background_3D = zeros(Nx, Ny, Nz);
% slice_nr = 14;
% background_3D(:,:,slice_nr) = mask_3D(:,:,slice_nr);
% ROI_img{r+5} = background_3D;
% I_roi{r+5} = find(background_3D);
% ROI_names{r+5} = ['Slice' num2str(slice_nr)];
% 
% N_ROI_REF = r+5;



% Get the boundaries of included mask voxels per slice, this is for
% visualisation purposes
bound = cell(1,N_ROIs);
for i = 1:N_ROIs
    bound{i} = cell(Nz,1);
    for k = 1:Nz
        bound{i}{k,1} = bwboundaries(squeeze(ROI_img{i}(:, :, k)));
    end    
end

