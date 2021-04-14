% fmrwhy_settings_template: Settings for the fmrwhy_workflow_qc pipeline

% Main data source: BIDS root folder
options.bids_dir = '/Users/jheunis/Desktop/ds002748';

% Subjects to run
options.subjects_output = 'all';

% Set template T1w image (set to '' if a single T1w image was collected)
options.anat_template_session = '';

% Set template for functional realignment purposes (if not needed, set to '')
options.template_task = 'rest';
options.template_session = '';
options.template_run = '';
options.template_echo = '';

% Sequence parameters
options.TR = '';
options.N_slices = 25;
options.Ndummies = 0;
options.Nscans = 100;
options.TE = []; % assume for all functional runs

% Settings for structFunc processing

% Settings for anatLocaliser processing
options.map_rois = 0;
% options.roi_orig_dir = '/Volumes/Stephan_WD/NEUFEPME_data_templates';
options.roi_orig_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_templates';
options.roi = struct;
% IMPORTANT: structure has to be named using the task name as in options.tasks: options.roi.(task).orig_fn
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

% Settings for basicFunc processing
options.fwhm = 5;
options.basicfunc_full = false; % if true, preprocessing will include all combinations of slice time correction, realignment and smoothing, useful for later analyses; if false, only include steps necessary for QC
options.include_stc = false;

% Settings for generateMultRegr routine
options.confounds.include_volterra = 1;
options.confounds.include_fd = 1;
options.confounds.include_tissue = 1;
options.confounds.include_physio = 0;

% generateMultRegr: framewise displacement
options.r = 50; % mm
options.FD_threshold = 0; % set as 0 to calculate with both standard thresholds 0.2 and 0.5 mm.

% generateMultRegr: PhysIO
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

% Settings for QC
options.theplot.intensity_scale = [-6 6];
options.theplot.include_physio = 0;
options.qc_overwrite_tissuecontours = false;
options.qc_overwrite_ROIcontours = false;
options.qc_overwrite_theplot = true;
options.qc_overwrite_statsoutput = false;

% Settings for first level analysis: steps to include/exclude
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
