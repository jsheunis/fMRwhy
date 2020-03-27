% -----------
% -----------
% TO FOCUS ON
% -----------
% -----------

% 1. QC pipeline:
%   - include task paradigm plots
% 2. Functional localisation pipeline:
% 3. VERY IMPORTANT INVESTIGATE USE OF DICM2NII EVERYWHERE WHEN LOADING NIFTI (AS OPPOSED TO SPMREADVOLS), DOES THIS INFLUENCE RESULTS?
%   -
%   -





% -----------------
% -----------------
% ISSUES TO ADDRESS
% -----------------
% -----------------

% --------------------------------------------------
% 1. ORIENTATION WHEN PLOTTING/loading BRAIN  SLICES
% --------------------------------------------------
% IMPORTANT: NEED TO FIGURE THIS OUT

% Info from https://brainder.org/2012/09/23/the-nifti-file-format/
% The most visible improvement of the nifti format over the previous analyze format is the ability
% to unambiguously store information orientation. The file standard assumes that the voxel coordinates
% refer to the center of each voxel, rather than at any of its corners. The world coordinate system is
% assumed to be ras: +x is Right, +y is Anterior and +z is Superior, which is precisely different than
% the coordinate system used in analyze, which is las. The format provides three different methods to map
% the voxel coordinates (i,j,k) to the world coordinates (x,y,z). The first method exists only to allow
% compatibility with the analyze format. The other two methods may coexist, and convey different coordinate
% systems. These systems are specified in the fields short qform_code and short sform_code, which can
% assume the values specified in the table:

% RAS+ = +Right +Anterior +Superior = +X +Y +Z

% Links:
% Useful description of everything: https://nipy.org/nibabel/nifti_images.html
% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=SPM;eb275f02.1908:
%   - says spm.mat reflects sform rows, which spm gives precedence to. Look at spm.private.mat0 for qform.


% Important note 21 March 2020:
% nii = nifti(func_fn) gives the same info as func_spm.private; e.g. nii.hdr = func_spm.private.hdr.

% Issue where qform and sform codes change after creating template 3d image from 4d timeseries:
% https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1808&L=SPM&P=R42327&K=2&X=E96B8262A335059C28&Y=j.s.heunis%40tue.nl


% func_spm.private:
% NIFTI object: 1-by-1
%            dat: [64×64×34 file_array]
%            mat: [4×4 double]
%     mat_intent: 'Aligned'
%           mat0: [4×4 double]
%    mat0_intent: 'Aligned'
%         timing: [1×1 struct]
%        descrip: 'Template functional volume'
% Stephan interpretation: mat0 = qform rows; intent = code.


% IMPORTANT: Decision on 23 March 2020: use code from https://github.com/xiangruili/dicm2nii in fmrwhy_util_readNifti:
% This reorients the nifti image in voxel space such that it is RAS+. For imagesc, this image then has to be
% rotated with 90 degrees, because for some reason (as of yet unknown to me) imagesc swops the x and y axes;
% or interprets these axes differently when provided in a matrix format.

% IMPORTANT: For now, only use nii_tools to load image data *for viewing*. Can expand later to viewing and processing and saving once I understand it better.

% IMPORTANT: NEED TO CHECK WHEREVER spm_read_vols and spm_vol ARE USED for viewing AND REPLACE IF REQUIRED!!!!!!!


% ---------------------------------------
% 2. HOW DEPENDENCIES ARE HANDLED
% ---------------------------------------

% e.g. currently have dependencies folder with dicm2nii, while physio is not included


% ---------------------------------------
% 3. HTML REPORTS
% ---------------------------------------

% - download bootstrap stuff to have them locally
% - embed images in html base64
% -




% -------------------
% -------------------
% EDUCATIONAL CONTENT
% -------------------
% -------------------

% System of linear equations: Y = Xb + e
%   - Y = [Ntime x Nvox] matrix =   timeseries of all brain voxels
%   - X = [Ntime x Nregr] matrix =  design matrix of regressors specifying our data model
%   - b = [Nregr x 1] matrix =      beta values corresponding to model regressors, to be solved
%   - e = [Ntime x Nvox] matrix =   error

% Solve system of linear equations:
% 1.  Y = Xb_hat
% ==> Multiply each side by transpose of X
% 2.  X'Y = (X'X)b_hat
% ==> If X has full column rank, X'X is invertible. Multiply each side by inverse of X'X
% 3.  inv(X'X)X'Y = inv(X'X)(X'X)b_hat
% 4.  inv(X'X)X'Y = b_hat
% 5. b_hat = estimate of b = same as least squares estimate
% In Matlab:
%Y = Xb + e
%b_hat = pinv(X)*Y
%Y_hat = Y - X*b_hat






