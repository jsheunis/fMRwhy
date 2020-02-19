% A custom workflow that does structFunc and basicFunc preprocessing and QC fof a single subject in the NEUFEP study

% Code steps:
% 1. Define template/default variables, directories and filenames
% 2. Specify


%--------------------------------------------------------------------------

% Load/create required defaults
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
bids_dir = '/Volumes/Stephan_WD/KempCourse_data_organised/bids';
sub = '017';
ses = '';
template_task = 'motor';
template_run = '1';
template_echo = '2';
task = 'motor';
run = '1';
echo = '2';

% BIDS structure values
% BIDS = spm_BIDS(bids_dir);
% subjects = spm_BIDS(BIDS,'subjects');
% sessions = spm_BIDS(BIDS,'sessions');
% runs = spm_BIDS(BIDS,'runs');
% tasks = spm_BIDS(BIDS,'tasks');
% types = spm_BIDS(BIDS,'types');
% modalities = spm_BIDS(BIDS,'modalities');

% Directory and content setup
deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
qc_dir = fullfile(deriv_dir, 'fmrwhy-qc');
sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);
sub_dir_qc = fullfile(qc_dir, ['sub-' sub]);
stats_dir = fullfile(deriv_dir, 'fmrwhy-stats');
sub_dir_stats = fullfile(stats_dir, ['sub-' sub]);
if ~exist(sub_dir_preproc, 'dir')
    mkdir(sub_dir_preproc)
    sub_dir_BIDS = fullfile(bids_dir, ['sub-' sub]);
    copyfile(sub_dir_BIDS, sub_dir_preproc)
end
if ~exist(sub_dir_qc, 'dir')
    mkdir(sub_dir_qc)
end
if ~exist(sub_dir_stats, 'dir')
    mkdir(sub_dir_stats)
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

%%
% -------
% STEP 0 -- Create functional template
% -------
tic;
% Define template filename
template_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_space-individual_bold.nii']);
% Create, if it does not exist
if ~exist(template_fn, 'file')
    disp(['Template funcional image does not exist yet. Creating now: ' template_fn]);
    functional0_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_echo-' template_echo '_bold.nii,1']);
    fmrwhy_util_saveNifti(template_fn, spm_read_vols(spm_vol(functional0_fn)), functional0_fn, 'Template functional volume', 0)
else
    disp(['Template funcional image exists: ' template_fn]);
end
toc;
%%
% -------
% STEP 1 -- Structural-functional preprocessing: fmrwhy_preproc_structFunc.m
% -------
struct_func_out_fns = [{coregest_anatomical_fn} probseg_fns transform_fns rall_fns mask_fns];
run_structFunc = 0;
for i = 1:numel(struct_func_out_fns)
    if ~exist(struct_func_out_fns{i}, 'file')
        disp(['Structural-funcional preprocessing output file does not exist yet: ' struct_func_out_fns{i}]);
        run_structFunc = 1;
    end
end
if run_structFunc
    disp('Running complete structural-funcional preprocessing pipeline')
    fmrwhy_preproc_structFunc(bids_dir, sub, ses, template_task, template_run, template_echo);
    disp('Complete!')
    disp('---')
else
    disp('Structural-funcional preprocessing already completed.')
    disp('---')
end
toc;
%%
% -------
% STEP 2 -- Basic functional preprocessing: fmrwhy_preproc_basicFunc.m
% -------
basic_func_out_fns = {motion_fn, afunctional_fn, rfunctional_fn, rafunctional_fn, sfunctional_fn, srfunctional_fn, srafunctional_fn, confounds_fn};
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
toc;
%%
% -------
% STEP 3 -- Quality control pipeline: fmrwhy_qc_.m
% -------
%qc_out_fns; % some file that is generated by the qc procedure (subject level or run level?)
run_qc = 1;
%for i = 1:numel(qc_out_fns)
%    if ~exist(qc_out_fns{i}, 'file')
%        disp(['Basic funcional preprocessing output file does not exist yet: ' qc_out_fns{i}]);
%        run_qc = 1;
%    end
%end

if run_qc
    disp('Running quality control pipeline.')
    fmrwhy_qc_run(bids_dir, sub, ses, template_task, template_run, template_echo);
    disp('Complete!')
    disp('---')
else
    disp('Basic quality control pipeline already completed.')
    disp('---')
end
toc;
%%
% -------
% STEP 4 -- 1st level analysis
% -------
func_dir_stats = fullfile(sub_dir_stats, 'func');
if ~exist(func_dir_stats, 'dir')
    mkdir(func_dir_stats)
end

% Set up statistical design parameters, based on task data
% Load multiple confound regressors
confounds_struct = tdfread(confounds_fn);
confounds_mat = struct2array(confounds_struct);
new_confounds_mat = confounds_mat(:,1:(end-1)); % ignore global signal regressor
if isempty(find(confounds_struct.framewise_displacement_censor))
    new_confounds_mat(:,8) = []; % remove censoring regressor if there are no censored volumes
end
new_confounds_mat(:,7) = []; % remove fd as a regressor
new_confounds_fn = fullfile(func_dir_stats, ['sub-' sub '_task-' task '_run-' run '_desc-GLM_regressors.txt']);

dlmwrite(new_confounds_fn, new_confounds_mat, 'delimiter', '\t', 'precision', '%1.7e')

sess_params = struct;
sess_params.timing_units = 'scans';
sess_params.timing_RT = 2;
sess_params.cond_name = 'Fingertapping';
sess_params.cond_onset = [11; 31; 51; 71; 91; 111; 131; 151; 171; 191];
sess_params.cond_duration = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10];
spm_specify1stlevel_jsh(func_dir_stats, srafunctional_fn, new_confounds_fn, sess_params)
%%
load([func_dir_stats filesep 'SPM.mat']);
%% ESTIMATE MODEL
fmrwhy_batch_estimate1stlevel(func_dir_stats)
%% SETUP TASK CONTRAST
[Ntt, Nregr] = size(SPM.xX.X);
contrast_params = struct;
contrast_params.weights = zeros(1, Nregr);
contrast_params.weights(1) = 1;
contrast_params.name = 'Fingertapping';
fmrwhy_batch_contrast1stlevel(func_dir_stats, contrast_params)
%% RUN RESULTS
fmrwhy_batch_threshold1stlevel(func_dir_stats)
[SPM,xSPM] = spm_getSPM(fullfile(func_dir_stats, 'SPM.mat'));
toc;