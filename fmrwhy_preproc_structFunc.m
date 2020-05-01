function fmrwhy_preproc_structFunc(bids_dir, sub, ses, task, run, echo, options)
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


% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% Update workflow params with subject functional derivative filenames
options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);




% -------
% STEP 1 -- Coregister (estimate) structural image to template functional image
% -------
if ~exist(options.coregest_anatomical_fn, 'file')
    disp('Coregistering T1w image to functional template')
    fmrwhy_batch_coregEst(options.anatomical_fn, options.template_fn, options.coregest_anatomical_fn);
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
for i = 1:numel(options.probseg_fns)
    if ~exist(options.probseg_fns{i}, 'file')
        disp(['Segmentation file does not exist yet: ' options.probseg_fns{i}]);
        run_seg = 1;
    end
end
if run_seg
    disp('Segmenting the coregistered T1w image into tissue compartments')
    fmrwhy_batch_segment(options.coregest_anatomical_fn, options.spm_dir, options.probseg_fns, options.transform_fns);
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
for i = 1:numel(options.rall_fns)
    if ~exist(options.rall_fns{i}, 'file')
        disp(['Resliced file does not exist yet: ' options.rall_fns{i}]);
        run_resl = 1;
    end
end
if run_resl
    disp('Resampling the coregistered T1w image and tissue compartments into template functional space')
    reslice_fns = {};
    reslice_fns{1} = options.coregest_anatomical_fn;
    for i = 2:7
        reslice_fns{i} = options.probseg_fns{i-1};
    end
    fmrwhy_batch_coregResl(reslice_fns, options.template_fn, options.rall_fns)
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
for i = 1:numel(options.mask_fns)
    if ~exist(options.mask_fns{i}, 'file')
        disp(['Mask does not exist yet: ' options.mask_fns{i}]);
        run_masks = 1;
    end
end
if run_masks
    disp('Computing and saving structural mask images in functional subject space')
    % Get binary 3D images for each tissue type, based on a comparison of
    % the probability value for each tissue type per voxel (after applying
    % a treshold on the probability values). Also combined GM, WM and CSF to get brain mask.
    [GM_img_bin, WM_img_bin, CSF_img_bin, brain_img_bin] = fmrwhy_util_createBinaryMasks(options.rgm_probseg_fn, options.rwm_probseg_fn, options.rcsf_probseg_fn, 0.5);
    % save masks to file: rtme_util_saveNifti(template_fn, img, new_fn, descrip)
    fmrwhy_util_saveNifti(options.gm_mask_fn, double(GM_img_bin), options.template_fn)
    fmrwhy_util_saveNifti(options.wm_mask_fn, double(WM_img_bin), options.template_fn)
    fmrwhy_util_saveNifti(options.csf_mask_fn, double(CSF_img_bin), options.template_fn)
    fmrwhy_util_saveNifti(options.brain_mask_fn, double(brain_img_bin), options.template_fn)
    disp('Complete!')
    disp('---')
else
    disp('Structural mask images in functional subject space already exist')
    disp('---')
end