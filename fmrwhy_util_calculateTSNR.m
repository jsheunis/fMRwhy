function output = fmrwhy_util_calculateTSNR(functional_fn, mask_fn, saveAs_fn, template_fn)
% Function to calculate tSNR from (masked) fMRI timeseries data

% INPUT:

% INPUT:

% OUTPUT:
output = struct;

vols = spm_vol(functional_fn);

nii = nii_tool('load', functional_fn);
data_4D = nii.img; % [Ni x Nj x Nk x Nt]
[Ni, Nj, Nk, Nt] = size(data_4D);
data_2D = reshape(data_4D, Ni*Nj*Nk, Nt); %[voxels, time]
data_2D = data_2D'; %[time, voxels]

% Then get the mask if fn supplied
if mask_fn ~= 0
    nii_mask = nii_tool('load', mask_fn);
    mask_3D = nii_mask.img; % [Ni x Nj x Nk]
    mask_2D = reshape(mask_3D, Ni*Nj*Nk, 1); %[voxels, 1]
    I_mask = find(mask_2D); %[voxels, 1]
    I_mask = I_mask'; %[1, voxels]
end

% Remove linear and quadratic trends from data, per voxel
data_2D_detrended = fmrwhy_util_detrend(data_2D, 2); %[time, voxels]
output.data_2D_detrended = data_2D_detrended;

% Calculate mean (ignore nan values)
data_2D_mean = nanmean(data_2D_detrended); %[1, voxels]
data_3D_mean = reshape(data_2D_mean, Ni, Nj, Nk);
output.data_3D_mean = data_3D_mean;

% Calculate standard deviation
data_2D_stddev = std(data_2D_detrended); %[1, voxels]
data_3D_stddev = reshape(data_2D_stddev, Ni, Nj, Nk);
output.data_3D_stddev = data_3D_stddev;

% Calculate tSNR
data_2D_tsnr = data_2D_mean./data_2D_stddev; %[1, voxels]
data_3D_tsnr = reshape(data_2D_tsnr, Ni, Nj, Nk);
output.data_3D_tsnr = data_3D_tsnr;

% Mask if required
if mask_fn ~= 0
    data_2D_tsnr_masked = zeros(size(data_2D_tsnr)); %[1, voxels]
    data_2D_tsnr_masked(I_mask) = tSNR_2D(I_mask); %[1, voxels]
    data_3D_tsnr_masked = reshape(data_2D_tsnr_masked, Ni, Nj, Nk);
    output.data_3D_tsnr_masked = data_3D_tsnr_masked;
end

% Save to file, if required
if saveAs_fn ~= 0
    if mask_fn ~= 0
        fmrwhy_util_saveNifti(saveAs_fn, data_3D_tsnr_masked, template_fn)
    else
        fmrwhy_util_saveNifti(saveAs_fn, data_3D_tsnr, template_fn)
    end
end
