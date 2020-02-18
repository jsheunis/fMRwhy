% STEPS TO RUN FOR PREPROCESSING AND QC
% ALSO SHOWS CORRESPONDING OUTPUT FILENAMES


% fmrwhy_preproc_structFunc(bids_dir, sub, ses, task, run, echo, defaults)
%       -








% Pre-processing: anatomical (RUN 1)
%
% 1. Anatomical to functional space coregistration, use middle echo first volume as template - SPM12 coregister estimate
% 2. Segment coregistered anatomical image into tissue components - SPM12 unified segmentation
%     1. Saves inverse transform from subject functional to MNI space
% 3. Coregister relevant regions of interest (from atlases in MNI space) to subject functional space using inverse transfromations
% 4. Reslice all to functional space grid (SPM reslice)
%
%
% Pre-processing: peripheral data (RUN 1)
%
% 1. Generate RETROICOR regressors from cardiac and respiratory traces of both runs (run 2 data to be used later) - PhysIO + Matlab
%
%
% Pre-processing: functional (RUN 1)
%
% 1. Task region localisation (using only middle echo [TE=28ms] timeseries):
%     1. Slice time correction
%     2. 3D volume realignment
%     3. Calculate framewise displacement from realignment params, select outliers using FD threshold (*which value or percentage?*)
%     4. Gaussian kernel smoothing (2*voxel size?)
%     5. GLM analysis incl:
%         1. AR(1) autoregressive filtering
%         2. Drift removal / high-pass (SPM cosine basis set)
%         3. Realignment params [+expansion?]
%         4. RETROICOR (+HRV, RTV?)
%         5. FD outlier binary regressor
%         6. *(global or tissue compartment signals???)*
%     6. Select t-stat peak within anatomically bound mask (from anatomy toolbox ROI)
%     7. Select N amount of voxels neighbouring peak voxel ==> ROI for real-time use
%
% 2. T2*, S0, tSNR calculation from `run1_BOLD_rest` dataset (*is this sensible, as opposed to using RUN 1 task data?*):
%     1. Slice time correction on all three echo timeseries
%     2. 3D volume realignment on middle echo timeseries
%     3. Apply rigid body transformations from middle echo realignment parameters to echo 1 and echo 3 timeseries
%     6. T2* and S0 estimation (*check steps of tedana*):
%         1. *How to mask?*
%         2. Calculate timeseries average
%         3. Estimate T2* and S0 using log-linear fit of mono-exponential decay model
%         4. *Threshold?*
%     4. *Drift removal?*
%     5. tSNR calculation:
%         1. *How to mask?*
%         2. Mean / stddev