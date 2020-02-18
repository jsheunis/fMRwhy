function stats = fmrwhy_qc_calculateTSNR(functional_fn, mask_fn)
% Function to calculate tSNR from (masked) fMRI timeseries data
%
% INPUT:
% MP                 - movement parameter matrix (Nt x 6)

% OUTPUT:

% First access and reshape the functional data: 4D to 2D
vols = spm_vol(functional_fn);
[N_i, N_j, N_k] = vols(1).dim;
N_t = numel(vols);
data_4D = spm_read_vols(vols);
data_2D = reshape(data_4D, N_i*N_j*N_k, N_t); %[voxels, time]
data_2D = data_2D'; %[time, voxels]

% Then get the mask
mask_vol = spm_vol(mask_fn);
mask_3D = spm_read_vols(mask_vol);
mask_2D = reshape(mask_3D, N_i*N_j*N_k, 1); %[voxels, 1]
I_mask = find(data_mask); %[voxels, 1]
I_mask = I_mask'; %[1, voxels]

% Remove linear trend per voxel
data_2D_detrended = fmrwhy_util_detrend(data_2D, 1); %[time, voxels]

% Calculate mean
data_2D_mean = nanmean(data_2D_detrended); %[1, voxels]

% Calculate standard deviation
data_2D_stddev = std(data_2D_detrended); %[1, voxels]

% Calculate tSNR
tSNR_2D = data_2D_mean./data_2D_stddev;
tSNR_brain = mean(tSNR_2D(I_mask));
tSNR_GM = mean(tSNR_2D(I_GM));
tSNR_WM = mean(tSNR_2D(I_WM));
tSNR_CSF = mean(tSNR_2D(I_CSF));