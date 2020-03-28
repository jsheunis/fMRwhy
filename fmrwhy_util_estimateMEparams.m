function output = fmrwhy_util_estimateMEparams(me_fns, TE, mask_fn, template_fn, saveAs_t2star_fn, saveAs_s0_fn)
% Function to estimate blabla

% INPUT:
%   me_fns      -
%   TE          - 1xN vector with echo times in ms

%

% OUTPUT:

% -------
% STEP 1: Check dimensions of parameters
% -------
Nme = numel(me_fns);
Nte = numel(TE);
if Nme ~= Nte
    disp(['Number of provided filenames (' num2str(Nme) ') do not match the number of provided echo times (' num2str(Nte)]);
    return;
end

% -------
% STEP 2: Initialise/set variables
% -------
output = struct;
me_data = struct;
fn = me_fns{1};
vols = spm_vol(fn);
Ni = vols(1).dim(1);
Nj = vols(1).dim(2);
Nk = vols(1).dim(3);
Nt = numel(vols);
mask_3D = spm_read_vols(spm_vol(mask_fn));
mask_2D = reshape(mask_3D, Ni*Nj*Nk, 1); %[voxels, 1]
I_mask = find(mask_2D); %[voxels, 1]
I_mask = I_mask'; %[1, voxels]
% TODO: implement adaptive mask similar to tedana? See: https://github.com/ME-ICA/tedana/blob/master/tedana/utils.py

% -------
% STEP 3: For each echo timeseries, read data, reshape data, detrend data, and calculate timeseries mean
% -------
for e = 1:Nme
    % Read functional timeseries data for echo = e
    data_4D = spm_read_vols(spm_vol(me_fns{e}));
    % Reshape to 2D matrix of time x voxels
    data_2D = reshape(data_4D, Ni*Nj*Nk, Nt); %[voxels, time]
    me_data(e).data_2D = data_2D'; %[time, voxels]
    % Remove linear and quadratic trends from data, per voxel
    me_data(e).data_2D_detrended = fmrwhy_util_detrend(me_data(e).data_2D, 2); %[time, voxels]
    % Calculate timeseries mean per voxel (ignoring NaN)
    me_data(e).data_mean = nanmean(me_data(e).data_2D_detrended); %[1, voxels]
end

% -------
% STEP 4: Estimate T2star and S0 parameters using log-linear regression of the mono-exponential decay equation
% -------
% Initialise parameter maps to zeros
T2star = zeros(Ni*Nj*Nk, 1); %[voxels, 1]
T2star_thresholded = T2star;
S0 = zeros(Ni*Nj*Nk, 1); %[voxels, 1]
S0_thresholded = S0;
% Create "design matrix" X
X = horzcat(ones(Nme, 1), -TE(:));
% Create "data matrix" Y
Y = [];
for e = 1:Nme
    Y = [Y; me_data(e).data_mean(I_mask)];
end
Y_beforemax = Y;
% Natural logarithm of zero or negative value is not defined
Y = max(Y,1e-11);
% Estimate "beta matrix" by solving set of linear equations
beta_hat = pinv(X)*log(Y);
% Calculate S0 and T2star from beta estimation
S0(I_mask) = exp(beta_hat(1, :));
T2star(I_mask) = 1./beta_hat(2, :);
% Now threshold the T2star (and S0?) values based on (criteria to be defined:
% TODO: perhaps derive from https://github.com/ME-ICA/tedana/blob/master/tedana/workflows/t2smap.py)
T2star_thresh_max = 500; % arbitrarily chosen, same as tedana. TODO: inspect to see if this yields sensible results
T2star_thresh_min = 0; % arbitrarily chosen, same as tedana
T2star_thresholded(I_mask) = T2star(I_mask); % values outside of the mask will stay zero
I_T2star_min = (T2star_thresholded(:) < T2star_thresh_min) ; % vector of voxels where T2star value is negative
T2star_thresholded(I_T2star_min) = 0; % if values inside mask are zero or negative, set them to threshold_min value
I_T2star_max =  T2star_thresholded(:) >= T2star_thresh_max; % vector of voxels in mask where T2star value is higher than threshold_max
T2star_thresholded(I_T2star_max) = 0; % set to zero, different from tedana where they set it to threshold_max value
S0_thresholded(I_mask) = S0(I_mask);
I_S0_nan = isnan(S0_thresholded(:));
S0_thresholded(I_S0_nan) = 0;
% Convert the estimated and corrected parameters to 3D matrices
T2star_3D = reshape(T2star, Ni, Nj, Nk);
T2star_3D_thresholded = reshape(T2star_thresholded, Ni, Nj, Nk);
S0_3D = reshape(S0, Ni, Nj, Nk);
S0_3D_thresholded = reshape(S0_thresholded, Ni, Nj, Nk);

% Save to file, if required
fmrwhy_util_saveNifti(saveAs_t2star_fn, T2star_3D_thresholded, template_fn, 'T2star map', 1);
fmrwhy_util_saveNifti(saveAs_s0_fn, S0_3D_thresholded, template_fn, 'S0 map', 0);

% Save output variables
output.T2star_3D = T2star_3D;
output.T2star_3D_thresholded = T2star_3D_thresholded;
output.S0_3D = S0_3D;
output.S0_3D_thresholded = S0_3D_thresholded;
output.me_data = me_data;
output.Y = Y;
output.X = X;
output.Y_beforemax = Y_beforemax;
output.beta_hat = beta_hat;
output.I_mask = I_mask;
output.I_T2star_min = I_T2star_min;
output.I_T2star_max = I_T2star_max;
output.I_S0_nan = I_S0_nan;



% ---------------------------------------------------------------
% Formula derivation for log-linear estimation of T2star and SO:
% ---------------------------------------------------------------
% S = S0*exp(-t/T2star) = S0*exp(-t*R2star)
% log(S) = log(S0*exp(-t/T2star)) = log(S0) - t*R2star
% ---------------------------------------------------------------
% [log(S(TE1))]        = [log(S0) - TE1*R2star]
% [log(S(TE2))]          [log(S0) - TE2*R2star]
% [log(S(TE3))]          [log(S0) - TE3*R2star]
%                     = [1 - TE1]     [log(S0)]
%                       [1 - TE2]  *  [R2star ]
%                       [1 - TE3]
% ---------------------------------------------------------------
% Y = X*b (+ e) ==> [log(S(TE))] = [1 - TE]  *  [log(S0)]
%                                               [R2star ]
% ---------------------------------------------------------------
% b_hat = pinv(X)*Y ==> [log(S0)] = pinv([1 - TE]) * [log(S(TE))]
%                       [R2star ]
% ---------------------------------------------------------------