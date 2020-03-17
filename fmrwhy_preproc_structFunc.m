function fmrwhy_preproc_structFunc(bids_dir, sub, ses, task, run, echo, opts)
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

% Setup fmrwhy bids directories on workflow level
fmrwhy_defaults_setupDerivDirs(bids_dir);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
fmrwhy_defaults_setupSubDirs(bids_dir, sub);

% Update workflow params with subject anatomical derivative filenames
opts = fmrwhy_defaults_subAnat(bids_dir, sub, opts);

% Update workflow params with subject functional derivative filenames
opts = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, opts);


% -------
% STEP 1 -- Coregister (estimate) structural image to template functional image
% -------
if ~exist(opts.coregest_anatomical_fn, 'file')
    disp('Coregistering T1w image to functional template')
    fmrwhy_batch_coregEst(opts.anatomical_fn, opts.template_fn, opts.coregest_anatomical_fn);
    disp('Complete!')
    disp('---')
else
    disp('T1w image already coregistered to functional template')
    disp('---')
end


% -------
% STEP 2 -- Segmentation of coregistered anatomical image into GM, WM, CSF, etc
% -------
run_seg = 0;
for i = 1:numel(opts.probseg_fns)
    if ~exist(opts.probseg_fns{i}, 'file')
        disp(['Segmentation file does not exist yet: ' opts.probseg_fns{i}]);
        run_seg = 1;
    end
end
if run_seg
    disp('Segmenting the coregistered T1w image into tissue compartments')
    fmrwhy_batch_segment(opts.coregest_anatomical_fn, opts.spm_dir, opts.probseg_fns, opts.transform_fns);
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
run_resl = 0;
for i = 1:numel(opts.rall_fns)
    if ~exist(opts.rall_fns{i}, 'file')
        disp(['Resliced file does not exist yet: ' opts.rall_fns{i}]);
        run_resl = 1;
    end
end
if run_resl
    disp('Resampling the coregistered T1w image and tissue compartments into template functional space')
    reslice_fns = {};
    reslice_fns{1} = opts.coregest_anatomical_fn;
    for i = 2:7
        reslice_fns{i} = opts.probseg_fns{i-1};
    end
    fmrwhy_batch_coregResl(reslice_fns, opts.template_fn, opts.rall_fns)
    disp('Complete!')
    disp('---')
else
    disp('Resampling of coregistered T1w image and tissue compartments previously completed.')
    disp('---')
end

% -------
% STEP 4 -- Construct GM, WM, CSF and whole brain (GM+WM+CSF) masks
% ------
run_masks = 0;
for i = 1:numel(opts.mask_fns)
    if ~exist(opts.mask_fns{i}, 'file')
        disp(['Mask does not exist yet: ' opts.mask_fns{i}]);
        run_masks = 1;
    end
end
if run_masks
    disp('Computing and saving structural mask images in functional subject space')
    % Get binary 3D images for each tissue type, based on a comparison of
    % the probability value for each tissue type per voxel (after applying
    % a treshold on the probability values). Also combined GM, WM and CSF to get brain mask.
    [GM_img_bin, WM_img_bin, CSF_img_bin, brain_img_bin] = fmrwhy_util_createBinaryMasks(opts.rgm_probseg_fn, opts.rwm_probseg_fn, opts.rcsf_probseg_fn, 0.5);
    % save masks to file: rtme_util_saveNifti(template_fn, img, new_fn, descrip)
    fmrwhy_util_saveNifti(opts.gm_mask_fn, GM_img_bin, opts.template_fn, 'GM mask', 1)
    fmrwhy_util_saveNifti(opts.wm_mask_fn, WM_img_bin, opts.template_fn, 'WM mask', 1)
    fmrwhy_util_saveNifti(opts.csf_mask_fn, CSF_img_bin, opts.template_fn, 'CSF mask', 1)
    fmrwhy_util_saveNifti(opts.brain_mask_fn, brain_img_bin, opts.template_fn, 'Brain mask', 1)
    disp('Complete!')
    disp('---')
else
    disp('Structural mask images in functional subject space already exist')
    disp('---')
end