function output = detrend4D(functional4D_fn)
% Function to detrend 4D fMRI data
%
% INPUT:
% functional4D_fn       - filename for 4D fMRI timeseries data (.nii)
% 
% OUTPUT: 
% F_4D_detrended        - 
%__________________________________________________________________________
% Copyright (C) Stephan Heunis 2018


% Get file info and convert 4D to 2D

V = spm_vol(functional4D_fn);
V1_3D = spm_read_vols(V(1));
sizeV = size(V);
[Ni, Nj, Nk] = size(V1_3D);
Nt = sizeV(1);

[x,y,z] = ind2sub(size(V1_3D),(1:Ni*Nj*Nk)');
dataT = spm_get_data(V,[x y z]');
F_2D = dataT';
F_4D = reshape(F_2D, Ni, Nj, Nk, Nt);

% Setup design matrix to include demeaned 1st, 2nd and 3rd order
% polynomial and single mean regressors
X_design = [ (1:Nt)' ((1:Nt).^2/(Nt^2))' ((1:Nt).^3/(Nt^3))'];
X_design = X_design - mean(X_design);
X_design = [X_design ones(Nt,1)];

% GLM to get beta coefficients
betas = X_design\F_2D';

% Detrend data
F_2D_detrended = F_2D' - X_design(:, 1:(end-1))*betas(1:(end-1), :);
F_2D_detrended = F_2D_detrended';

% 4D matrix
F4D_detrended = reshape(F_2D_detrended, Ni, Nj, Nk, Nt);

% Output
output.F_2D_detrended = F_2D_detrended;
output.F4D_detrended = F4D_detrended;

