function output = fmrwhy_qc_createMaskMontages(bids_dir, sub, savefigs)
% Function to create montages of tissue mask contours overlaid on
% assumes standard fmrwhy-preproc directory structure
% assumes fmrwhy_preproc_structFunc.m has been run successfully


% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, []);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% Get functional template
template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
template_anat_fn = options.rcoregest_anatomical_fn;
[p, frm, rg, dim] = fmrwhy_util_readOrientNifti(template_anat_fn);
template_img = p.nii.img;

% Get anatomical masks in individual functional space
masks = fmrwhy_util_loadOrientMasks(bids_dir, sub);
mask_images = {masks.GM_mask_3D, masks.WM_mask_3D, masks.CSF_mask_3D, masks.brain_mask_3D};
mask_names = {'Grey matter', 'White matter', 'Cerebrospinal fluid', 'Whole brain'};
mask_montage_fns = {'_GM_mask_montage', '_WM_mask_montage', '_CSF_mask_montage', '_brain_mask_montage'};
for i = 1:numel(mask_montage_fns)
    mask_montage_fns{i} = fullfile(options.anat_dir_qc, ['sub-' sub mask_montage_fns{i} '.png']);
end

% Structure to save output
output = struct;

% Overlay montages
for i = 1:numel(mask_images)
    mask_img = mask_images{i};
    if savefigs
        fmrwhy_util_createOverlayMontage(template_img, mask_img, 9, 1, '', 'gray', 'off', 'max', [], [255,0,0], mask_montage_fns{i});
    else
        fmrwhy_util_createOverlayMontage(template_img, mask_img, 9, 1, '', 'gray', 'off', 'max', [], [255,0,0], 0);
    end
end