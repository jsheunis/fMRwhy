function fmrwhy_preproc_structFunc(bids_dir, sub, ses, task, run, echo)
%--------------------------------------------------------------------------

% Copyright statement....

%--------------------------------------------------------------------------
% DEFINITION
%--------------------------------------------------------------------------

% Function for pre-real-time anatomical/structural to functionap
% preprocessing for a single subject. Steps include coregistering
% structural image to initial functional image,segmenting the coregistered
% structural image into tissue types, and reslicing the segments to the
% functional resolution image grid. Makes use of spm12 batch routines.
% If spm12 batch parameters are not explicitly set, defaults are assumed.

% STEPS:
% 1. Anatomical to functional space coregistration, use middle echo first volume rest run 1 as template - SPM12 coregister estimate
% 2. Segment coregistered anatomical image into tissue components - SPM12 unified segmentation
%     - Saves inverse transform from subject functional to MNI space
% 3. Reslice all to functional space grid (SPM reslice)
% 4. Create tissue compartment and whole brain masks

% INPUT:
% functional0_fn     - filename of initial pre-real-time 3D functional volume template
% anatomical_fn     - filename of T1-weighted structural volume
% spm_dir           - SPM12 directory

% OUTPUT:

%--------------------------------------------------------------------------

disp('---')
disp('*** Running fmrwhy_preproc_structFunc ***')
disp('---')
disp('---')

% Load/create required defaults
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
template_task = 'motor'; % changed for fingertapping experiment. TODO: change back. and update functioning.
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

% Get anatomical file
anatomical_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_T1w.nii']);

% Get functional template
template_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_space-individual_bold.nii']);

% Initialise output filenames
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
% Outputs after reslicing segments and coregistered T1w image
rcoregest_anatomical_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-coregEstResl_T1w.nii']);
rgm_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rGM_probseg.nii']);
rwm_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rWM_probseg.nii']);
rcsf_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rCSF_probseg.nii']);
rbone_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rbone_probseg.nii']);
rsoft_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rsoft_probseg.nii']);
rair_probseg_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-rair_probseg.nii']);
% Outputs after creating masks
gm_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-GM_mask.nii']);
wm_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-WM_mask.nii']);
csf_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-CSF_mask.nii']);
brain_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-brain_mask.nii']);
% Structure to save outputs
output = struct;

% -------
% STEP 1 -- Coregister (estimate) structural image to template functional image
% -------
if ~exist(coregest_anatomical_fn, 'file')
    disp('Coregistering T1w image to functional template')
    fmrwhy_batch_coregEst(anatomical_fn, template_fn, coregest_anatomical_fn);
    disp('Complete!')
    disp('---')
else
    disp('T1w image already coregistered to functional template')
    disp('---')
end


% -------
% STEP 2 -- Segmentation of coregistered anatomical image into GM, WM, CSF, etc
% -------
probseg_fns = {gm_probseg_fn, wm_probseg_fn, csf_probseg_fn, bone_probseg_fn, soft_probseg_fn, air_probseg_fn};
transform_fns = {indiv_to_mni_fn, mni_to_indiv_fn};
run_seg = 0;
for i = 1:numel(probseg_fns)
    if ~exist(probseg_fns{i}, 'file')
        disp(['Segmentation file does not exist yet: ' probseg_fns{i}]);
        run_seg = 1;
    end
end
if run_seg
    disp('Segmenting the coregistered T1w image into tissue compartments')
    fmrwhy_batch_segment(coregest_anatomical_fn, spm_dir, probseg_fns, transform_fns);
    disp('Complete!')
    disp('---')
else
    disp('Tissue compartment images already exist (i.e. segmentation previously completed)')
    disp('---')
end

% -------
% STEP 3 -- Reslice all of the above to functional-resolution image grid
% TODO: check if it is necessary to copy files to temporary filenames before reslicing, or not
% -------
rprobseg_fns = {rgm_probseg_fn, rwm_probseg_fn, rcsf_probseg_fn, rbone_probseg_fn, rsoft_probseg_fn, rair_probseg_fn};
rall_fns = {rcoregest_anatomical_fn, rgm_probseg_fn, rwm_probseg_fn, rcsf_probseg_fn, rbone_probseg_fn, rsoft_probseg_fn, rair_probseg_fn};
run_resl = 0;
for i = 1:numel(rall_fns)
    if ~exist(rall_fns{i}, 'file')
        disp(['Resliced file does not exist yet: ' rall_fns{i}]);
        run_resl = 1;
    end
end
if run_resl
    disp('Resampling the coregistered T1w image and tissue compartments into template functional space')
    reslice_fns = {};
    reslice_fns{1} = coregest_anatomical_fn;
    for i = 2:7
        reslice_fns{i} = probseg_fns{i-1};
    end
    fmrwhy_batch_coregResl(reslice_fns, template_fn, rall_fns)
    disp('Complete!')
    disp('---')
else
    disp('Resampling of coregistered T1w image and tissue compartments previously completed.')
    disp('---')
end

% -------
% STEP 4 -- Construct GM, WM, CSF and whole brain (GM+WM+CSF) masks
% ------
mask_fns = {gm_mask_fn, wm_mask_fn, csf_mask_fn, brain_mask_fn};
run_masks = 0;
for i = 1:numel(mask_fns)
    if ~exist(mask_fns{i}, 'file')
        disp(['Mask does not exist yet: ' mask_fns{i}]);
        run_masks = 1;
    end
end
if run_masks
    disp('Computing and saving structural mask images in functional subject space')
    % Get binary 3D images for each tissue type, based on a comparison of
    % the probability value for each tissue type per voxel (after applying
    % a treshold on the probability values). Also combined GM, WM and CSF to get brain mask.
    [GM_img_bin, WM_img_bin, CSF_img_bin, brain_img_bin] = fmrwhy_util_createBinaryMasks(rgm_probseg_fn, rwm_probseg_fn, rcsf_probseg_fn, 0.5);
    % save masks to file: rtme_util_saveNifti(template_fn, img, new_fn, descrip)
    fmrwhy_util_saveNifti(gm_mask_fn, GM_img_bin, template_fn, 'GM mask', 1)
    fmrwhy_util_saveNifti(wm_mask_fn, WM_img_bin, template_fn, 'WM mask', 1)
    fmrwhy_util_saveNifti(csf_mask_fn, CSF_img_bin, template_fn, 'CSF mask', 1)
    fmrwhy_util_saveNifti(brain_mask_fn, brain_img_bin, template_fn, 'Brain mask', 1)
    disp('Complete!')
    disp('---')
else
    disp('Structural mask images in functional subject space already exist')
    disp('---')
end