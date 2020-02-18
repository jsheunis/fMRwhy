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
