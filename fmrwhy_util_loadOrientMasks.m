function masks = fmrwhy_util_loadOrientMasks(bids_dir, sub)
% Load data for masks that were already calculated;
% assumes standard fmrwhy-preproc directory structure
% assumes fmrwhy_preproc_structFunc.m has been run successfully

% This function loads mask data (using fmrwhy_util_readNifti derived from dicm2nii) FOR PLOTTING purposes.
% It includes reorientation for RAS+ display.
% If pure calculation purposes are intended, rather use fmrwhy_util_loadMasks!

deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);

masks = struct;

masks.GM_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-GM_mask.nii']);
masks.WM_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-WM_mask.nii']);
masks.CSF_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-CSF_mask.nii']);
masks.brain_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-brain_mask.nii']);

[p, frm, rg, dim] = fmrwhy_util_readOrientNifti(masks.GM_mask_fn);
masks.GM_mask_3D = p.nii.img; % [Ni x Nj x Nk]
[p, frm, rg, dim] = fmrwhy_util_readOrientNifti(masks.WM_mask_fn);
masks.WM_mask_3D = p.nii.img; % [Ni x Nj x Nk]
[p, frm, rg, dim] = fmrwhy_util_readOrientNifti(masks.CSF_mask_fn);
masks.CSF_mask_3D = p.nii.img; % [Ni x Nj x Nk]
[p, frm, rg, dim] = fmrwhy_util_readOrientNifti(masks.brain_mask_fn);
masks.brain_mask_3D = p.nii.img; % [Ni x Nj x Nk]

[Ni, Nj, Nk] = size(p.nii.img);

masks.GM_mask_2D = reshape(masks.GM_mask_3D, Ni*Nj*Nk, 1); % [Ni*Nj*Nk x 1]
masks.WM_mask_2D = reshape(masks.WM_mask_3D, Ni*Nj*Nk, 1); % [Ni*Nj*Nk x 1]
masks.CSF_mask_2D = reshape(masks.CSF_mask_3D, Ni*Nj*Nk, 1); % [Ni*Nj*Nk x 1]
masks.brain_mask_2D = reshape(masks.brain_mask_3D, Ni*Nj*Nk, 1); % [Ni*Nj*Nk x 1]

masks.GM_mask_I = find(masks.GM_mask_2D); % [Ngmvoxels x 1]
masks.WM_mask_I = find(masks.WM_mask_2D); % [Nwmvoxels x 1]
masks.CSF_mask_I = find(masks.CSF_mask_2D); % [Ncsfvoxels x 1]
masks.brain_mask_I = find(masks.brain_mask_2D); % [Nbrainvoxels x 1]

masks.GM_mask_I = masks.GM_mask_I'; % [1 x Ngmvoxels]
masks.WM_mask_I = masks.WM_mask_I'; % [1 x Nwmvoxels]
masks.CSF_mask_I = masks.CSF_mask_I'; % [1 x Ncsfvoxels]
masks.brain_mask_I = masks.brain_mask_I'; % [1 x Nbrainvoxels]

masks.field_names = {'GM', 'WM', 'CSF', 'brain'};
