function masks = fmrwhy_util_loadMasks(bids_dir, sub)
% Load data for masks that were already calculated;
% assumes standard fmrwhy-preproc directory structure
% assumes fmrwhy_preproc_structFunc.m has been run successfully

deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);

masks = struct;

masks.GM_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-GM_mask.nii']);
masks.WM_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-WM_mask.nii']);
masks.CSF_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-CSF_mask.nii']);
masks.brain_mask_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_space-individual_desc-brain_mask.nii']);

masks.GM_mask_3D = spm_read_vols(spm_vol(masks.GM_mask_fn)); % [Ni x Nj x Nk]
masks.WM_mask_3D = spm_read_vols(spm_vol(masks.WM_mask_fn)); % [Ni x Nj x Nk]
masks.CSF_mask_3D = spm_read_vols(spm_vol(masks.CSF_mask_fn)); % [Ni x Nj x Nk]
masks.brain_mask_3D = spm_read_vols(spm_vol(masks.brain_mask_fn)); % [Ni x Nj x Nk]

vol = spm_vol(masks.brain_mask_fn);
Ni = vol.dim(1);
Nj = vol.dim(2);
Nk = vol.dim(3);

masks.GM_mask_2D = reshape(masks.GM_mask_3D, Ni*Nj*Nk, 1); % [Ni*Nj*Nk x 1]
masks.WM_mask_2D = reshape(masks.WM_mask_3D, Ni*Nj*Nk, 1); % [Ni*Nj*Nk x 1]
masks.CSF_mask_2D = reshape(masks.CSF_mask_3D, Ni*Nj*Nk, 1); % [Ni*Nj*Nk x 1]
masks.brain_mask_2D = reshape(masks.brain_mask_3D, Ni*Nj*Nk, 1); % [Ni*Nj*Nk x 1]

masks.GM_mask_I = find(masks.GM_mask_2D); % [Ni*Nj*Nk x 1]
masks.WM_mask_I = find(masks.WM_mask_2D); % [Ni*Nj*Nk x 1]
masks.CSF_mask_I = find(masks.CSF_mask_2D); % [Ni*Nj*Nk x 1]
masks.brain_mask_I = find(masks.brain_mask_2D); % [Ni*Nj*Nk x 1]

masks.GM_mask_I = masks.GM_mask_I'; % [1 x Ni*Nj*Nk]
masks.WM_mask_I = masks.WM_mask_I'; % [1 x Ni*Nj*Nk]
masks.CSF_mask_I = masks.CSF_mask_I'; % [1 x Ni*Nj*Nk]
masks.brain_mask_I = masks.brain_mask_I'; % [1 x Ni*Nj*Nk]

masks.field_names = {'GM', 'WM', 'CSF', 'brain'};
