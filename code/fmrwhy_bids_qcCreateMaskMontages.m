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
mask_montage_fns = {'_GM_mask_montage', '_WM_mask_montage', '_CSF_mask_montage', '_brain_mask_montage'};
if isempty(options.anat_template_session)
    [filename, filepath] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'ext', '_T1w.nii');
else
    [filename, filepath] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'ses', options.anat_template_session, 'ext', '_T1w.nii');
end
for i = 1:numel(mask_montage_fns)
    mask_montage_fns{i} = fullfile(options.qc_dir, filepath, filename);
    mask_montage_fns{i} = strrep(mask_montage_fns{i}, '_T1w.nii', mask_montage_fns{i});
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