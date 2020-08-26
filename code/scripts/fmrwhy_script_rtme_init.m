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

% TODO: NOTE - EVERYTHING IS READ IN AND SAVED WITH SPM_VOL (ETC) AND NOT USING NII_TOOL

% -------
% STEP 3: Directories, files, parameters (from fmrwhy output) to be used during rt processing
% -------
Ne = options.Ne; % number of echoes per volume
% Templates
functional0_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
funcref_spm = spm_vol(functional0_fn);
funcref_3D  = spm_read_vols(funcref_spm);
[Nx, Ny, Nz] = size(funcref_3D);
t2star_fn = fullfile(options.sub_dir_me, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_desc-MEparams_t2star.nii']);
t2star_img = spm_read_vols(spm_vol(t2star_fn));
s0_fn = fullfile(options.sub_dir_me, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_desc-MEparams_s0.nii']);
s0_img = spm_read_vols(spm_vol(s0_fn));
tsnr_data = zeros([Nx Ny Nz Ne]);
tsnr_fn = {};
for e = 1:Ne
    tsnr_fn{e} = fullfile(options.sub_dir_me, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_echo-' num2str(e) '_desc-rapreproc_tsnr.nii']);
    tsnr_data(:,:,:,e) = spm_read_vols(spm_vol(tsnr_fn{e}));
end
% Raw bold volumes to process in real-time
boldts_fn = {};
for e = 1:Ne
    boldts_fn{e} = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' num2str(e) '_bold.nii']);
    boldts_3Dvolumes = dir(fullfile(options.sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_echo-' num2str(e) '_bold*.nii']));
    if numel(boldts_3Dvolumes) ~= options.Nscans
        disp('Splitting 4D scan to 3D volumes...')
        Vo = spm_file_split(boldts_fn{e}, options.sub_dir_rt);
    end
end

% Get tissue and whole brain masks
masks = fmrwhy_util_loadMasksSPM(bids_dir, sub);
N_tissuemasks = 4;
% Get task-based regions of interest (resliced and in subject functional space)
N_taskROIs = 2; %numel(options.roi.(task).rroi_fn);
N_ROIs = N_taskROIs + N_tissuemasks;
ROI_fns = cell(1,N_ROIs);
ROI_img = cell(1,N_ROIs);
I_roi = cell(1,N_ROIs);
ROI_names = cell(1,N_ROIs);
noFWE_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-2_noFWEp001e20']);
if strcmp(task, 'motor')
    roi_anat_fn = fullfile(options.sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rleftMotor_roi.nii']);
    anat_desc = 'rleftMotor';
    if strcmp(run, '1')
        funcAnat_desc = 'FingerTappingOverlapsleftMotornoFWE';
    else
        funcAnat_desc = 'MentalFingerTappingOverlapsleftMotornoFWE';
    end
    roi_funcAnat_fn = fullfile(noFWE_dir_stats, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-' funcAnat_desc '_roi.nii']);

else
    roi_anat_fn = fullfile(options.sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rbilateralAmygdala_roi.nii']);
    anat_desc = 'rbilateralAmygdala';
    if strcmp(run, '1')
        funcAnat_desc = 'Faces>ShapesOverlapsbilateralAmygdalanoFWE';
    else
        funcAnat_desc = 'MentalEmotionOverlapsbilateralAmygdalanoFWE';
    end
    roi_funcAnat_fn = fullfile(noFWE_dir_stats, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-' funcAnat_desc '_roi.nii']);
end

% Anatomical ROI
r = 1;
ROI_fns{r} = roi_anat_fn;
ROI_img{r} = spm_read_vols(spm_vol(ROI_fns{r}));
ROI_img{r}(isnan(ROI_img{r})) = 0;
ROI_img{r} = (ROI_img{r} > 0.1);
I_roi{r} = find(ROI_img{r}(:)); % find(ROI_img{r}(:) & GM_img_bin(:)); TODO, investigate which masking to use here, if any
ROI_names{r} = anat_desc;
% Functionally localised and anatomically constrained ROI
r = 2;
ROI_fns{r} = roi_funcAnat_fn;
ROI_img{r} = spm_read_vols(spm_vol(ROI_fns{r}));
ROI_img{r}(isnan(ROI_img{r})) = 0;
ROI_img{r} = (ROI_img{r} > 0.1);
I_roi{r} = find(ROI_img{r}(:)); % find(ROI_img{r}(:) & GM_img_bin(:)); TODO, investigate which masking to use here, if any
ROI_names{r} = funcAnat_desc;
% Gray matter ROI
ROI_names{r+1} = 'GM';
ROI_img{r+1} = masks.GM_mask_3D;
I_roi{r+1} = masks.GM_mask_I;
% White matter ROI
ROI_names{r+2} = 'WM';
ROI_img{r+2} = masks.WM_mask_3D;
I_roi{r+2} = masks.WM_mask_I;
% CSF ROI
ROI_names{r+3} = 'CSF';
ROI_img{r+3} = masks.CSF_mask_3D;
I_roi{r+3} = masks.CSF_mask_I;
% Whole brain ROI
ROI_names{r+4} = 'brain';
ROI_img{r+4} = masks.brain_mask_3D;
I_roi{r+4} = masks.brain_mask_I;
% Main brain mask
I_mask = masks.brain_mask_I;
N_maskvox = numel(I_mask);


% -------
% STEP 3: Define params for real-time experiment
% -------
% Volume dimensions, and reference image
Nvox = Nx*Ny*Nz;
ROI_native = true; % are ROIs already in native space with matching resolution?
TE = options.TE; % Echo times in ms
use_echo = 0; % specify if only a single echo should be used for real-time; set to 0 to use all echoes
template_echo = 2; % echo to be used for motion parameters, for ME
T2star_thresh = 120; % threshold for maximum T2star after estimation (= T2* of CSF at 3T, Cesar)
voxel_size = [3.5 3.5 3.5];
smoothing_kernel    = [7 7 7];
Nt =   options.Nscans; % NrOfVolumes % VolumesNumber
N_skip = 0; % nrSkipVol
N_start = N_skip + 1;
Ndyn = Nt - N_skip;
TR = 2;
TR_skip = 2; % amount of TRs to skip at start of baseline block to exclude effect of HRF on baseline mean calculation
NF_cond = 2; % location of task/nf condition in design matrix (SPM structure)
lCond = 2; % number of conditions (should actually be derived form SPM mat)
% TODO: automate task conditions
task_onsets = [11; 31; 51; 71; 91; 111; 131; 151; 171; 191];
task_durations = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10];
baseline_onsets = [1; 21; 41; 61; 81; 101; 121; 141; 161; 181; 201];
baseline_durations = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10];
% For the visualisation purposes:
Nslice = 25;
rotateDir = [0 0 1];
rotateVal = 1;
rotateDeg = 270;
% Get the boundaries of included mask voxels per slice
bound = cell(1,N_ROIs);
for i = 1:N_ROIs
    bound{i} = cell(Nz,1);
    for k = 1:Nz
        bound{i}{k,1} = bwboundaries(squeeze(ROI_img{i}(:, :, k)));
    end
end


% -------
% STEP 4: Define parameters for real-time processing (from OpenNFT)
% -------
% Processing steps to include: AR(1), cGLM, iGLM
isIGLM = false;
iglmAR1 = false;
isRegrIGLM = false;
isMotionRegr = true;
isHighPass = true;
isLinRegr = true;
isPhysRegr = false;
cglmAR1 = false;

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

for sig = 1:7
    tmp_posMin{sig}(1:N_ROIs) = 0;
    tmp_posMax{sig}(1:N_ROIs) = 0;
end

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

% Protocol parameters (derived/calculated)
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
% Calculate HRF-convolved task and baseline designs. The HRF (hemodynamic
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
    basFct = arRegr_opennft(aAR1, basFct);
end
nrBasFct = size(basFct,2);

% Other OpenNFT values without explanation, but necessary
spmMaskTh = 80*ones(Ndyn,1)'; %80*ones(size(SPM.xM.TH));% mean(SPM.xM.TH)*ones(size(SPM.xM.TH)); % SPM.xM.TH;
pVal = .01;
statMap3D_iGLM = [];
tContr = [0; 1]; % Baseline; Task]  TODO

%% Data setup

% Get motion realignment template data and volume
dicomInfoVox   = sqrt(sum((funcref_spm.mat(1:3,1:3)).^2));
fwhm = smoothing_kernel ./ dicomInfoVox;
A0=[];x1=[];x2=[];x3=[];wt=[];deg=[];b=[];
% MULTI-ECHO TODO
R = struct;
R(1,1).mat = funcref_spm.mat;
R(1,1).dim = funcref_spm.dim;
R(1,1).Vol = funcref_3D;
MP = nan(Ndyn, 6);


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
    spmDesign = arRegr_opennft(aAR1, tmpSpmDesign);
end
mf = [];
npv = 0;
statMapCreated = 0;
nrRegrDesign = size(spmDesign,2);
% number of regressors of no interest to correct with cGLM
nrRegrToCorrect = 8; % 6 MC regressors, linear trend, constant

