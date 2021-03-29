function FD_measures = fmrwhy_qc_calculateFD(MP, r, FD_threshold)
% Function to calculate framewise displacement (FD) and related measures.
%
% INPUT:
% MP                 - movement parameter matrix (Nt x 6)
% r                  - radius to use for displacement derived from
%                      rotational movement parameters
% FD_threshold       - threshold (in mm) that defines outlier volumes
%                    - if threshold == 0, use both predefined trhesholds of 0.2 and 0.5 mm
%
% OUTPUT:
% FD_measures        - structure with filenames and data
%__________________________________________________________________________
% Copyright (C) Stephan Heunis 2018

% Define variables
FD_measures = struct;

% First demean and detrend (linear) the movement parameters
MP_corrected = fmrwhy_util_detrend(MP, 1);

% Then transform rotational parameters to linear translations (small angle assumption)
MP_mm = MP_corrected;
MP_mm(:,4:6) = MP_mm(:,4:6)*r; % 50mm from Power 2017; 80 mm from QAP

% Calculate FD and related measures
MP_diff = [zeros(1, 6); diff(MP_mm)];
FD_measures.FD = sum(abs(MP_diff),2);

if FD_threshold == 0
    FD_measures.FD_outliers_regr02 = FD_measures.FD>=0.2;
    FD_measures.FD_outliers_ind02 = find(FD_measures.FD_outliers_regr02);
    FD_measures.FD_outliers_regr05 = FD_measures.FD>=0.5;
    FD_measures.FD_outliers_ind05 = find(FD_measures.FD_outliers_regr05);
else
    FD_measures.FD_outliers_regr = FD_measures.FD>=FD_threshold;
    FD_measures.FD_outliers_ind = find(FD_measures.FD_outliers_regr);
end

FD_measures.FD_sum = sum(FD_measures.FD);
FD_measures.FD_mean = mean(FD_measures.FD);