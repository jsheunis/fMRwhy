function output = fmrwhy_bids_qcCreateMaskMontages(bids_dir, sub, savefigs, options)
% Function to create montages of tissue mask contours overlaid on
% assumes standard fmrwhy-preproc directory structure
% assumes fmrwhy_preproc_structFunc.m has been run successfully

%-------------
% Run script
%-------------

% Get template background image
template_anat_fn = options.rcoregest_anatomical_fn;
[p, frm, rg, dim] = fmrwhy_util_readOrientNifti(template_anat_fn);
template_img = p.nii.img;

% Get anatomical masks in individual functional space
masks = fmrwhy_util_loadOrientMasks(bids_dir, sub, options);
mask_images = {masks.GM_mask_3D, masks.WM_mask_3D, masks.CSF_mask_3D, masks.brain_mask_3D};
mask_names = {'Grey matter', 'White matter', 'Cerebrospinal fluid', 'Whole brain'};
mask_desc = {'GM', 'WM', 'CSF', 'brain'};

for i = 1:numel(mask_images)
    if isempty(options.anat_template_session)
        [filename, filepath] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', mask_desc{i}, 'ext', '_mask_montage.png');
    else
        [filename, filepath] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'ses', options.anat_template_session, 'space', 'individual', 'desc', mask_desc{i}, 'ext', '_mask_montage.png');
    end
    mask_montage_fns{i} = fullfile(options.qc_dir, filepath, filename);
end

% Structure to save output
output = struct;

% Overlay montages
for i = 1:numel(mask_images)
    mask_img = mask_images{i};
    if savefigs
        fmrwhy_util_createOverlayMontage(template_img, mask_img, 7, 1, '', 'gray', 'off', 'max', [], [255,0,0], mask_montage_fns{i});
    else
        fmrwhy_util_createOverlayMontage(template_img, mask_img, 7, 1, '', 'gray', 'off', 'max', [], [255,0,0], 0);
    end
end