% ------------------------
% fmrwhy_settings_template
% ------------------------

% ----------
% Section 01
% ----------

% Main data source: BIDS root folder
options.bids_dir = '';

% ----------
% Section 02
% ----------

% Subjects to run. If all ==> 'all'. If a subset, include the subject identifier as a cell array.
options.subjects_output = 'all'; % e.g.: options.subjects_output = {'001', '003', '021'};

% ----------
% Section 03
% ----------
% Realignment uses the SPM12 2step procedure: realignment to first, followed by realignment to mean
% You can select the level at which realignment occurs. Options:
% per_task:
% - all runs of a task are included (in order) in the realignment procedure
% per_run:
% - a single run is included in the realignment procedure; this is repeated for all runs.
% to_template:
% - a single run is included in the realignment procedure; it is realigned to a template specified below. this is repeated for all runs.
options.realignment_type = 'per_task'; % per_task / per_run / to_template
options.realignment_template_session = '';
options.realignment_template_task = '';
options.realignment_template_run = '';
options.realignment_template_echo = '';
% (if not needed, set to '')

% Set to which reference the coregistration should be done. Options:
% per_task (i.e. multiple coregistrations):
% - coregistration is done to each task
% - the functional template is either specified by options.template_run (mean image)
%   or it is taken as the mean image of all runs of the same task (and session)
% per_run (i.e. multiple coregistrations):
% - coregistration is done to each run
% - the functional template is taken as the mean image of the specific run (options.template_run is ignored)
% to_template (i.e. a single coregistration):
options.coreg_type = 'per_task'; % per_task / per_run / to_template
% Set the template for functional data, used for coregistration
options.coreg_template_session = '';
options.coreg_template_task = '';
options.coreg_template_run = '';
options.coreg_template_echo = '';
% (if not needed, set to '')

% Set the T1w image that will be used for anatomical-to-functional registration steps
% Default behaviour is that session-specific T1w images are used for coregistration of functional
% data in the same session (options.anat_template_session = '')
% Custom behaviour is introduced if a template session is selected (e.g. options.anat_template_session = '1')
% In the latter case the template session T1w image is copied to other sessions (and renamed accordingly)
% in order to simplify processing steps and avoid duplication+overwriting of results. This
% is then followed by the same procedure as session-specific coregistration.
% Copying is only done in the derivatives directory, BIDS dataset data are left untouched.
% This has obvious (possible bad) implications for interpreting the preprocessed data,
% but it is deemed fine if we remain aware of this.
options.anat_template_session = '';


% ----------
% Section 04
% ----------

% Define sequence parameters, this should be known to the user.
% fMRwhy does not yet support deriving these parameters from the BIDS json files
options.TR = 2;
options.N_slices = 34;
options.Ndummies = 5; % Specify number of dummies, even if they have already been excluded from the BOLD timeseries images. This is important for the TAPAS PhysIO steps when processing cardiac and respiratory data.
options.Nscans = 210; % The number of scans in the BOLD timeseries images, excluding dummies.
options.TE = [14 28 42]; % The echo times if multi-echo. fMRwhy assumes this is the same for all functional runs. If not multi-echo, set to [].

% ----------
% Section 05 - Settings for anatLocaliser processing
% ----------

% Should the QC pipeline also map ROIs to the subject space. Yes = 1; No = 0. If 0, the rest of this section can be ignored.
options.map_rois = 0;
% Specify directory with roi files, assumed to be in MNI152 space
options.roi_orig_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_templates';
options.roi = struct;
% Specify the roi filenames, names, descriptions
% IMPORTANT: structure has to be named using the task name as specified in the BIDS data structure: options.roi.(task).orig_fn
options.roi.motor.orig_fn = {fullfile(options.roi_orig_dir, 'Left_Motor_4a_4p.nii')
                             fullfile(options.roi_orig_dir, 'Right_Motor_4a_4p.nii')}; % Raw ROI filenames

options.roi.motor.name = {'Left Motor', 'Right Motor'}; % For plots and strings
options.roi.motor.desc = {'leftMotor', 'rightMotor'}; % For BIDS file naming (after normalisation  to functional space)

options.roi.emotion.orig_fn = {fullfile(options.roi_orig_dir, 'Bilateral_Amygdala_allregions.nii')
                               fullfile(options.roi_orig_dir, 'Left_Amygdala_allregions.nii')
                               fullfile(options.roi_orig_dir, 'Right_Amygdala_allregions.nii')}; % Raw ROI filenames

options.roi.emotion.name = {'Bilateral Amygdala', 'Left Amygdala', 'Right Amygdala'}; % For plots and strings
options.roi.emotion.desc = {'bilateralAmygdala', 'leftAmygdala', 'rightAmygdala'}; % For BIDS file naming (after normalisation  to functional space)

% options.roi.(task).roi_fn = ROIs in subject space (not resliced)
% options.roi.(task).rroi_fn = resliced ROIs in subject space

% ----------
% Section 06 - Settings for basicFunc processing
% ----------

options.fwhm = 7; % FWHM of the spatial smoothing kernel, typically twice the voxel size, in mm
options.basicfunc_full = false; % if true, preprocessing will include all combinations of slice time correction, realignment and smoothing, useful for later analyses; if false, only include steps necessary for QC
options.include_stc = false; % include slice timing correction?

% Settings for generateMultRegr routine; specifies which regressors to generate for GLM analysis
options.confounds.include_volterra = 1; % to generate the volterra expansion of the 6 realignment parameters ==> 1
options.confounds.include_fd = 1; % to generate framewise displacement ==> 1
options.confounds.include_tissue = 1; % to generate signals from GM, WM, CSF compartments ==> 1
options.confounds.include_physio = 1; % if cardiac and respiration data are available ==> 1

% generateMultRegr: framewise displacement. No need to change these defaults unless for very specific and well-motivated reasons.
options.r = 50; % mm
options.FD_threshold = 0; % set as 0 to calculate with both standard thresholds 0.2 and 0.5 mm.

% ----------
% Section 07 - Settings for physiology processing
% ----------

% Parameters required by TAPAS PhysIO
% For a better understanding of how to set these parameters, read the TAPAS PhysIO documentation
% Parameters to set: sampling_interval, align_scan, onset_slice, cardiac_modality.
options.physio.options.cardiac_fn = '';
options.physio.options.respiration_fn = '';
options.physio.options.vendor = 'BIDS';
options.physio.options.sampling_interval = 0.002; % 500 Hz ==> Philips wired acquisition
options.physio.options.align_scan = 'last';
options.physio.options.Nslices = options.N_slices;
options.physio.options.TR = options.TR; % in seconds
options.physio.options.Ndummies = options.Ndummies; % include, even if these are not included in the fMRI timeseries data exported from the scanner
options.physio.options.Nscans = options.Nscans;
options.physio.options.onset_slice = 1;
options.physio.options.cardiac_modality = 'PPU';
options.physio.options.output_multiple_regressors_fn = 'PhysIO_multiple_regressors.txt'; % text file name
options.physio.options.level = 0; % verbose.level = 0 ==> do not generate figure outputs
options.physio.options.fig_output_file = ''; % unnecessary if verbose.level = 0, but still initialized here

% ----------
% Section 08 - Settings for QC processing and reporting
% ----------

% Settings for QC
% No need to change these defaults unless for very specific and well-motivated reasons.
options.theplot.intensity_scale = [-6 6];
options.qc_overwrite_tissuecontours = true;
options.qc_overwrite_ROIcontours = true;
options.qc_overwrite_theplot = false;
options.qc_overwrite_statsoutput = true;
options.qc_report_runs = {'task-rest_run-1', 'task-fingerTapping', 'task-emotionProcessing', 'task-rest_run-2', 'task-fingerTappingImagined', 'task-emotionProcessingImagined'};
options.qc_dataset_name = 'rt-me-fMRI';
options.qc_anat_res = '1x1x1 mm (100x100x100 voxels)';
options.qc_func_res = '3.5x3.5x3.5 mm (64x64x34 voxels)';
options.qc_func_acq = 'Multi-echo (TE = 14,28,42 ms), SENSE = 2.5';
options.qc_func_runs = 'rest_run-1, fingerTapping, emotionProcessing, rest_run-2, fingerTappingImagined, emotionProcessingImagined';

% ONLY FOR TASK DATA:

% ----------
% Section 09 - Settings for first level analysis
% ----------

% Settings for first level analysis: steps to include/exclude
% No need to change these defaults unless for very specific and well-motivated reasons.
options.firstlevel.tmap_montages = true;
options.firstlevel.anat_func_roi = true;

% Settings for first level analysis: task-motor
options.firstlevel.motor.run1.sess_params.timing_units = 'secs';
options.firstlevel.motor.run1.sess_params.timing_RT = 2;
options.firstlevel.motor.run1.sess_params.cond_names = {'FingerTapping'};
options.firstlevel.motor.run2.sess_params.timing_units = 'secs';
options.firstlevel.motor.run2.sess_params.timing_RT = 2;
options.firstlevel.motor.run2.sess_params.cond_names = {'MentalFingerTapping'};

% Settings for first level analysis: task-emotion
options.firstlevel.emotion.run1.sess_params.timing_units = 'secs';
options.firstlevel.emotion.run1.sess_params.timing_RT = 2;
options.firstlevel.emotion.run1.sess_params.cond_names = {'Faces', 'Shapes'};
options.firstlevel.emotion.run2.sess_params.timing_units = 'secs';
options.firstlevel.emotion.run2.sess_params.timing_RT = 2;
options.firstlevel.emotion.run2.sess_params.cond_names = {'MentalEmotion'};

% Settings for plotting task conditions
onset = [11; 31; 51; 71; 91; 111; 131; 151; 171; 191];
duration = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10];
options.firstlevel.motor.run1.plot_params.cond_onset = onset;
options.firstlevel.motor.run1.plot_params.cond_duration = duration;
options.firstlevel.motor.run2.plot_params.cond_onset = onset;
options.firstlevel.motor.run2.plot_params.cond_duration = duration;
options.firstlevel.emotion.run2.plot_params.cond_onset = onset;
options.firstlevel.emotion.run2.plot_params.cond_duration = duration;
onset = [12; 32; 52; 72; 92; 112; 132; 152; 172; 192];
duration = [9; 9; 9; 9; 9; 9; 9; 9; 9; 9];
options.firstlevel.emotion.run1.plot_params.cond_onset = onset;
options.firstlevel.emotion.run1.plot_params.cond_duration = duration;

% Settings for first level analysis: glm regressors to include
options.firstlevel.glm_regressors.trans_rot = true;
options.firstlevel.glm_regressors.trans_rot_derivative1 = true;
options.firstlevel.glm_regressors.trans_rot_power2 = false;
options.firstlevel.glm_regressors.trans_rot_derivative1_power2 = false;
options.firstlevel.glm_regressors.framewise_displacement_censor02 = false;
options.firstlevel.glm_regressors.framewise_displacement_censor05 = false;
options.firstlevel.glm_regressors.dvars_censor = false; % not yet implemented
options.firstlevel.glm_regressors.std_dvars_censor = false; % not yet implemented
options.firstlevel.glm_regressors.grey_matter = false;
options.firstlevel.glm_regressors.white_matter = false;
options.firstlevel.glm_regressors.csf = true;
options.firstlevel.glm_regressors.global_signal = false;
% Order of included retroicor regressors; if 0 ==> exclude
options.firstlevel.glm_regressors.retroicor_c = 2; % cardiac, max 6
options.firstlevel.glm_regressors.retroicor_r = 2; % respiratory, max 8
options.firstlevel.glm_regressors.retroicor_cxr = 0; % interaction, max 4
options.firstlevel.glm_regressors.hrv = false;
options.firstlevel.glm_regressors.rvt = false;

% Settings for first level analysis: task-motor
options.firstlevel.motor.run1.contrast_params.consess{1}.tcon.name = 'FingerTapping';
options.firstlevel.motor.run1.contrast_params.consess{1}.tcon.weights = [1];
options.firstlevel.motor.run1.contrast_params.consess{1}.tcon.sessrep = 'none';
options.firstlevel.motor.run2.contrast_params.consess{1}.tcon.name = 'MentalFingerTapping';
options.firstlevel.motor.run2.contrast_params.consess{1}.tcon.weights = [1];
options.firstlevel.motor.run2.contrast_params.consess{1}.tcon.sessrep = 'none';

% Settings for first level analysis: task-emotion
options.firstlevel.emotion.run1.contrast_params.consess{1}.tcon.name = 'Faces';
options.firstlevel.emotion.run1.contrast_params.consess{1}.tcon.weights = [1 0];
options.firstlevel.emotion.run1.contrast_params.consess{1}.tcon.sessrep = 'none';
options.firstlevel.emotion.run1.contrast_params.consess{2}.tcon.name = 'Shapes';
options.firstlevel.emotion.run1.contrast_params.consess{2}.tcon.weights = [0 1];
options.firstlevel.emotion.run1.contrast_params.consess{2}.tcon.sessrep = 'none';
options.firstlevel.emotion.run1.contrast_params.consess{3}.tcon.name = 'Faces>Shapes';
options.firstlevel.emotion.run1.contrast_params.consess{3}.tcon.weights = [1 -1];
options.firstlevel.emotion.run1.contrast_params.consess{3}.tcon.sessrep = 'none';
options.firstlevel.emotion.run2.contrast_params.consess{1}.tcon.name = 'MentalEmotion';
options.firstlevel.emotion.run2.contrast_params.consess{1}.tcon.weights = [1];
options.firstlevel.emotion.run2.contrast_params.consess{1}.tcon.sessrep = 'none';

% matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Patients > Control';
% matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [-1 1];
% matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
