function ROI_time_series = fmrwhy_util_getROItimeseries(functional_fn, mask_fn, roi_fn)
%
% ------------------------------------------------------------------------ %
% functional_fn = 4D timeseries;
% mask_fn = binary brain mask; exclude voxels outside of this mask
% roi_fn = filename of ROI
% ------------------------------------------------------------------------ %

% ------ %
% STEP 1 %
% ------ %
% Load timeseries data, get image and time dimensions, transform data
nii = nii_tool('load', functional_fn);
data_4D = double(nii.img);
[Ni, Nj, Nk, Nt] = size(data_4D); % [Ni x Nj x Nk x Nt]
data_2D = reshape(data_4D, Ni*Nj*Nk, Nt); %[voxels, time]
data_2D = data_2D'; %[time, voxels]

% Load mask data
mask_nii = nii_tool('load', mask_fn);
mask_3D = mask_nii.img; % [Ni x Nj x Nk]
mask_2D = reshape(mask_3D, Ni*Nj*Nk, 1); % [Ni*Nj*Nk x 1]
mask_I = find(mask_2D); % [Nmaskvoxels x 1]
mask_I = mask_I'; % [1 x Nmaskvoxels]

% Load ROI data
roi_nii = nii_tool('load', roi_fn);
roi_3D = roi_nii.img; % [Ni x Nj x Nk]
roi_3D = fmrwhy_util_createBinaryImg(roi_3D, 0.1);
roi_2D = reshape(roi_3D, Ni*Nj*Nk, 1); % [Ni*Nj*Nk x 1]
roi_I = find(roi_2D); % [Nroivoxels x 1]
roi_I = roi_I'; % [1 x Nroivoxels]

% Combine ROI and mask
roi_masked_2D = roi_2D & mask_2D;
roi_masked_I = find(roi_masked_2D); % [Nroivoxels x 1]
roi_masked_I = roi_masked_I'; % [1 x Nroivoxels]

% Remove linear and quadratic trend per voxel
data_2D_detrended = fmrwhy_util_detrend(data_2D, 2); %[time, voxels]

% Calculate mean
data_2D_mean = nanmean(data_2D_detrended); %[1, voxels]

% Demean
data_2D_demeaned = data_2D_detrended - data_2D_mean; %[time, voxels]

% Calculate standard deviation
data_2D_stddev = std(data_2D_detrended); %[1, voxels]

% Calculate percentage signal change: [I(t) - mean(I)]/mean(I)*100
data_2D_psc = 100*(data_2D_detrended./data_2D_mean) - 100;
data_2D_psc(isnan(data_2D_psc)) = 0;
F_2D_psc = data_2D_psc';

% Order voxels
ROI_signals = F_2D_psc(roi_masked_I, :);
ROI_time_series = mean(ROI_signals, 1);