% --------------------------------- %
% --------------------------------- %
% - NEUFEPME WORKFLOW DESCRIPTION - %
% --------------------------------- %
% --------------------------------- %

% ----------------------------------------------------------------
% STEP -42: We hold these truths to be completely non-self-evident
% ----------------------------------------------------------------

% 1) rest-run-1, in particular, volume 1 echo 2, is used as the main functional template for various reasons:
% - Necessary to calculate baseline estimates of T2*, S0, tSNR and other values (we don'' use task-based due to additional variance)
% - Need a motion realignment and ROI templates before real-time is started (i.e. for rest2, motor2, emotion2)
% - Most similar to single echo acquisition, i.e. allows comparison of single-echo and multi-echo pipelines and outcomes
% 2) Motion realignment of all multi-echo runs: realign 2nd echo timeseries to template, apply same rigid body transformations to other echoe timeseries
% 3)
% 4)
% 5)
% 6)
% 7)

% -----------------------------------
% STEP 0: Data
% -----------------------------------
% FUNCTION: describe data
% -----------------------------------

% 1. T1w
% 2. Rest1
% 3. Motor1
% 4. Emotion1
% 5. Rest2
% 6. Motor2
% 7. Emotion2


% ------------------------ %
% ------------------------ %
% STEP 0: Data preparation %
% ------------------------ %
% ------------------------ %

% 1. Define template functional run, volume and echo


% ----------------------------------- %
% ----------------------------------- %
% STEP 1: Pre-processing: anatomical  %
% ----------------------------------- %
% ----------------------------------- %
% FUNCTION: fmrwhy_preproc_structFunc %
% ----------------------------------- %
% ----------------------------------- %

% -----------
% DESCRIPTION
% -----------
% This workflow completes multiple steps that are necessary to bring the subject's anatomical image,
% tissue segmentations, and MNI-space regions of interest into the exact space and grid of the subject's
% functional runs. A template functional image has to be provided.

% -----
% STEPS
% -----
% 1. Anatomical to functional space coregistration to template functional volume - SPM12 coregister estimate
% 2. Segment coregistered anatomical image into tissue components - SPM12 unified segmentation
%   - Saves formward and inverse transforms from subject functional to MNI space
% 3. Construct GM, WM, CSF and whole brain (GM+WM+CSF) masks
% 4. Coregister relevant regions of interest (from Anatomy toolbox atlases in MNI space) to subject functional space using inverse transfromations
% 5. Reslice all to functional space grid (SPM reslice)

% ------------------
% COMMENTS and TODOs
% ------------------
% - add ROI functionality
%   - Left_Amygdala_allregions.nii
%   - Right_Amygdala_allregions.nii
%   - Left_Motor_4a_4p.nii
%   - Right_Motor_4a_4p.nii
% - consider splitting off ROI functionality into separate workflow step (i.e. not included in basicStruct)


% ---------------------------------- %
% ---------------------------------- %
% STEP 2: Pre-processing: functional %
% ---------------------------------- %
% ---------------------------------- %
% FUNCTION: fmrwhy_preproc_basicFunc %
% ---------------------------------- %
% ---------------------------------- %

% -----------
% DESCRIPTION
% -----------
% This workflow executes basic functional preprocessing on all functional runs in the subject's dataset.
% It serves two main purposes: preprocessing data for (1) the quality control (QC) workflow, and (2) the
% task analysis (GLM) workflow.

% -----
% STEPS
% -----

% 1. Estimate 3D volume realignment parameters from raw data of middle (template) echo (only estimate, no reslicing yet)
% 2. Slice timing correction on raw data (==> 'afunctional')
% 3. 3D volume realignment:
%   - of raw data (==> 'rfunctional')
%   - of slice timing corrected data (==> 'rafunctional')
% 4. Spatial smoothing:
%   - of raw data (==> 'sfunctional')
%   - of slice timing corrected data (==> 'safunctional')
%   - of realigned data (==> 'srfunctional')
%   - of slice timing corrected and realigned data (==> 'srafunctional')
% 5. Generate multiple regressors for GLM analysis and QC (all outputs to individual tsv files, and one combined tsv file)
%   - 3D realignment parameters (using run timeseries and template volume)
%   - framewise displacement (FD) and FD censoring (using realignment parameters and threshold defaults)
%   - tissue compartment signals (from coregistered anatomical segmentations)
%   - RETROICOR and HRV+RVT (from respiratory and ppu signals, using PhysIO)

% ------------------
% COMMENTS and TODOs
% ------------------
% - automate the handling of multi-echo processing
%   - automatically detect multi-echo and run necessary preproc steps on extra echoes
%   - which echo timeseries to use as default for realignment parameters
% - automate the workflow for BIDS compatibility, i.e. how/when to preprocess sessions, runs, etc
% - consider splitting off multiple regressor functionality into separate workflow step (i.e. not included in basicFunc)
% - for multiple regressors: add DVARS, and whatever else seems useful
% - for multiple regressors: change workflow such that all regressors are calculated, and then study-specific workflow/defaults has to specify which ones to use:
%   - realignment params: full Volterra expansion (squares and derivatives)
%   - FD and censoring: calculate using different radii and censor thresholds
%   - PhysIO: use cardiac (3rd order, 6 terms), respiratory (4th order, 8 terms), interaction (4 terms), hrv and rvt regressors (as current)



% ----------------------------------- %
% ----------------------------------- %
% STEP 3: QC: anatomical + functional %
% ----------------------------------- %
% ----------------------------------- %
% FUNCTION: fmrwhy_qc_run             %
% ----------------------------------- %
% ----------------------------------- %

% -----------
% DESCRIPTION
% -----------
% This workflow executes quality control processes on anatomical and functional data.
% For anatomical quality control, montages of the tissue segmentation masks+outlines overlaid on top of the
% functional template volume are generated. These help to see if coregistration of the T1w image was sucessful, and if
% segmentation was successful. For functional quality control, steps include generating montages of a run timeseries
% tSNR image, mean timeseries image, standard deviation image and more, which all help to visually spot artefacts or
% signal shortcomings. Functional QC also generates timeseries images, mainly "the plot" with FD, tissue compartment,
% physiological and other signal traces.

% -----
% STEPS
% -----

% ------------------
% COMMENTS and TODOs
% ------------------
% - How to define exclusion criteria based on QC outputs? FD? DVARS
% - Stil need to incorporate DVARS, stdDVARS, GCOR, (what else)?
% - Which timeseries should be used for calculating statistical measures like tSNR, etc?



% -------------------------------------- %
% -------------------------------------- %
% STEP 4: Task localisation - functional %
% -------------------------------------- %
% -------------------------------------- %
% FUNCTION: ?                         %
% -------------------------------------- %
% -------------------------------------- %

% -----------
% DESCRIPTION
% -----------
% Task region localisation using middle echo timeseries

% -----
% STEPS
% -----
% 0. First set up the following in workflow defaults file:
%   - 1st level design per task (if this needs to be calculated from BIDS events files, do that here)
%   - contrasts per task (if this needs to be calculated somehow, do that here)
% 1. Access output from basicFunc processing, specifically 'srafunctional', to be used in GLM analysis
% 2. From generateMultRegr output, extract regressors required for this particular workflow and write to tsv and text file:
%   - Regressors: Realignment params [only parameters and derivatives, i.e. 12 traces]
%   - Regressors: FD outlier binary regressor (**threshold to be decided, Cesar="for extreme motion")
%   - Regressors: Tissue compartment signals (only CSF)
%   - Regressors: RETROICOR (cardiac x 2, respiratory x 2, interaction x 0) and HRV, RVT
% 3. Specify 1st level design for GLM analysis (fmrwhy_batch_specify1stlevel), include:
%   - AR(1) autoregressive filtering (standard)
%   - Drift removal / high-pass (SPM cosine basis set, standard)
%   - Multiple regressors text file (from above)
% 4. Estimate the model (fmrwhy_batch_estimate1stlevel)
% 5. Specify and run contrasts (fmrwhy_batch_contrast1stlevel)
% 6. Threshold the results (fmrwhy_batch_threshold1stlevel)
% 7. Save xSPM.mat file which enables accessing the thresholded clusters
% 8. Thresholding and selection criteria???
%   - N amount of voxels neighbouring peak voxel

% ------------------
% COMMENTS and TODOs
% ------------------
% - need to set up a standardised and intuitive way to distinguish and use app-level defaults and settings for a workflow
% - need an intuitive way to automate generation of task design and contrast from defaults to GLM
% - need an intuitive way of including/excluding nuisance regressors in task-specific analysis workflow



% ---------------------------------------- %
% ---------------------------------------- %
% STEP x: Combine localisations and select %
% ---------------------------------------- %
% ---------------------------------------- %
% FUNCTION: ?                              %
% ---------------------------------------- %
% ---------------------------------------- %

%   - Select t-stat peak within anatomically bound mask (from anatomy toolbox ROI)
%   - Select N amount of voxels neighbouring peak voxel ==> ROI for real-time use




% ------------------------------------------ %
% ------------------------------------------ %
% STEP x: Multi-echo quantitative processing %
% ------------------------------------------ %
% ------------------------------------------ %
% FUNCTION: ?                                %
% ------------------------------------------ %
% ------------------------------------------ %

% 2. T2*, S0, tSNR calculation from `run1_BOLD_rest` dataset
%     1. Slice time correction on all three echo timeseries
%     2. 3D volume realignment on middle echo timeseries
%     3. Apply rigid body transformations from middle echo realignment parameters to echo 1 and echo 3 timeseries
%     4. T2* and S0 estimation (*check steps of tedana*):
%         1. *How to mask?*
%         2. Calculate timeseries average
%         3. Estimate T2* and S0 using log-linear fit of mono-exponential decay model
%         4. *Threshold?*
%     5. *Drift removal?*
%     6. tSNR calculation:
%         1. *How to mask?*
%         2. Which timeseries to use?
%         3. Mean / stddev



% ------------------------------------------ %
% ------------------------------------------ %
% IMPORTANT PROBLEMS THAT STILL NEED SOLVING %
% ------------------------------------------ %
% ------------------------------------------ %

% 1.    In quality report, montages of timeseries average images saved to file, then read from file and then plotted
%       seem to look different than montage of template volume. Investigate if saving nifti changes voxel positions.
%       This important to check when using masks and precalculated quantitative maps and comparing with real-time data.