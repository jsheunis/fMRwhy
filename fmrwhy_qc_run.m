function fmrwhy_qc_run(bids_dir, sub, ses, task, run)

% Software dependences:
% - Matlab vX
% - SPM12 r7771
% - (TAPAS PhysIO Toolbox vX)

% Data required in order for this function to run as intended:
% - T1w (raw data)
% - Single run of fMRI timeseries (raw data)
% - Head motion/movement parameters derived from unprocessed data (i.e. from realignment of raw fMRI timeseries)
% - (T1w coregistered to template functional volume (from fmrwhy_preproc_structFunc.m))
% - (Segmentations of GM, WM, CSF in template functional volume space (from fmrwhy_preproc_structFunc.m))

% Code steps:
% 1. Get BIDS directory of run
% 2. Create relevant derivative directories (qc and preproc) if it doesn't exist yet
% 3. QC steps for fmri runs in subject folder:
%   - Head motion/movement parameters derived from unprocessed data (i.e. from realignment of raw fMRI timeseries)
%   - Framewise displacement
%   - Statistical measures / images (tsnr, variance, std, psc, DVARS)
%   - The plot
%   - PhysIO quality metrics and plots

% INPUT
% - flag to generate pictures or not
% - flag to generate nii images or not


% Load/create required defaults
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
template_task = 'rest';
template_run = '1';
template_echo = '2';

% Directory and content setup
deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
qc_dir = fullfile(deriv_dir, 'fmrwhy-qc');
sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);
sub_dir_qc = fullfile(qc_dir, ['sub-' sub]);
func_dir_qc = fullfile(sub_dir_qc, 'func');
if ~exist(func_dir_qc, 'dir')
    mkdir(func_dir_qc)
end

% Grab anatomical image
anatomical_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_T1w.nii']);

% Grab functional timeseries filename,
functional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_bold.nii']);

% Grab template filename
template_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_space-individual_bold.nii']);


% ------------------------
% SECTION A: ANATOMICAL QC
% ------------------------

% -------
% STEP 1: Contours of tissue masks on mean EPI / template EPI (/ anatomical image?)
% -------
montage_data = fmrwhy_qc_createMaskMontages(bids_dir, sub, 1)
% -------
% STEP 2: Contours of anatomical ROIs on mean EPI / template EPI (/ anatomical image?)
% -------


% ------------------------
% SECTION B: FUNCTIONAL QC
% ------------------------

% -------
% STEP 1: Grab multiple regressors
% -------
motion_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_motion.tsv']);
framewise_displacement_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_fd.tsv']);
tissue_regr_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_tissue.tsv']);
physio_regr_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_physio.tsv']);
confounds_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_regressors.tsv']);

% -------
% STEP 2: Calculate and generate statistical measures / images (tsnr, variance, std, psc, DVARS)
% -------
stats = fmrwhy_qc_calculateStats(bids_dir, sub, ses, task, run)

% -------
% STEP 3: Plot The Plot
% -------
output = fmrwhy_qc_createThePlot(bids_dir, sub, ses, task, run)