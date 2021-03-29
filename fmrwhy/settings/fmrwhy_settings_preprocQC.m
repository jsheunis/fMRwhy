function settings = fmrwhy_settings_preprocQC(bids_dir, settings)

% Settings that are specific to a workflow to process a dataset

% BIDS structure values
%settings.BIDS = spm_BIDS(bids_dir);
%settings.subjects = spm_BIDS(BIDS,'subjects');
%settings.sessions = spm_BIDS(BIDS,'sessions');
%settings.runs = spm_BIDS(BIDS,'runs');
%settings.tasks = spm_BIDS(BIDS,'tasks');
%settings.types = spm_BIDS(BIDS,'types');
%settings.modalities = spm_BIDS(BIDS,'modalities');

% Set template for functional realignment purposes
settings.template_task = 'rest';
settings.template_session = '';
settings.template_run = '1';
settings.template_echo = '2';

% Sequence parameters
settings.TR = 2;
settings.N_slices = 34;
settings.Ndummies = 5;
settings.Nscans = 210;
settings.TE = [14 28 42];
settings.Ne = numel(settings.TE);

% Dataset parameters
settings.sessions = {[]};
settings.Nsessions = numel(settings.sessions);
settings.tasks = {'rest', 'motor', 'emotion'};
settings.Ntasks = numel(settings.tasks);
settings.runs = {'1', '2'};
settings.Nruns = numel(settings.runs);

% Settings for structFunc processing

% Settings for anatLocaliser processing
settings.map_rois = 1;
%settings.roi_orig_dir = '/Volumes/Stephan_WD/NEUFEPME_data_templates';
settings.roi_orig_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_templates';
settings.roi = struct;
% IMPORTANT: structure has to be named using the task name as in settings.tasks: settings.roi.(task).orig_fn
settings.roi.motor.orig_fn = {fullfile(settings.roi_orig_dir, 'Left_Motor_4a_4p.nii'),
                                fullfile(settings.roi_orig_dir, 'Right_Motor_4a_4p.nii')}; % Raw ROI filenames

settings.roi.motor.name = {'Left Motor', 'Right Motor'}; % For plots and strings
settings.roi.motor.desc = {'leftMotor', 'rightMotor'}; % For BIDS file naming (after normalisation  to functional space)

settings.roi.emotion.orig_fn = {fullfile(settings.roi_orig_dir, 'Bilateral_Amygdala_allregions.nii'),
                                fullfile(settings.roi_orig_dir, 'Left_Amygdala_allregions.nii'),
                                fullfile(settings.roi_orig_dir, 'Right_Amygdala_allregions.nii')}; % Raw ROI filenames

settings.roi.emotion.name = {'Bilateral Amygdala', 'Left Amygdala', 'Right Amygdala'}; % For plots and strings
settings.roi.emotion.desc = {'bilateralAmygdala', 'leftAmygdala', 'rightAmygdala'}; % For BIDS file naming (after normalisation  to functional space)

%settings.roi.(task).roi_fn = ROIs in subject space (not resliced)
%settings.roi.(task).rroi_fn = resliced ROIs in subject space


% Settings for basicFunc processing
settings.fwhm = 7;

% Settings for generateMultRegr routine
settings.confounds.include_volterra = 1;
settings.confounds.include_fd = 1;
settings.confounds.include_tissue = 1;
settings.confounds.include_physio = 1;

% generateMultRegr: framewise displacement
settings.r = 50; % mm
settings.FD_threshold = 0; % set as 0 to calculate with both standard thresholds 0.2 and 0.5 mm.

% generateMultRegr: PhysIO
settings.physio.options.cardiac_fn = '';
settings.physio.options.respiration_fn = '';
settings.physio.options.vendor = 'BIDS';
settings.physio.options.sampling_interval = 0.002; % 500 Hz ==> Philips wired acquisition
settings.physio.options.align_scan = 'last';
settings.physio.options.Nslices = settings.N_slices;
settings.physio.options.TR = settings.TR; % in seconds
settings.physio.options.Ndummies = settings.Ndummies; % include, even if these are not included in the fMRI timeseries data exported from the scanner
settings.physio.options.Nscans = settings.Nscans;
settings.physio.options.onset_slice = 1;
settings.physio.options.cardiac_modality = 'PPU';
settings.physio.options.output_multiple_regressors_fn = 'PhysIO_multiple_regressors.txt'; % text file name
settings.physio.options.level = 0; % verbose.level = 0 ==> do not generate figure outputs
settings.physio.options.fig_output_file = ''; % unnecessary if verbose.level = 0, but still initialized here


% Settings for QC
settings.theplot.intensity_scale = [-6 6];
settings.qc_overwrite_tissuecontours = true;
settings.qc_overwrite_ROIcontours = true;
settings.qc_overwrite_theplot = false;
settings.qc_overwrite_statsoutput = true;


% Settings for first level analysis: steps to include/exclude
settings.firstlevel.tmap_montages = true;
settings.firstlevel.anat_func_roi = true;

% Settings for first level analysis: task-motor
settings.firstlevel.motor.run1.sess_params.timing_units = 'secs';
settings.firstlevel.motor.run1.sess_params.timing_RT = 2;
settings.firstlevel.motor.run1.sess_params.cond_names = {'FingerTapping'};
settings.firstlevel.motor.run2.sess_params.timing_units = 'secs';
settings.firstlevel.motor.run2.sess_params.timing_RT = 2;
settings.firstlevel.motor.run2.sess_params.cond_names = {'MentalFingerTapping'};

% Settings for first level analysis: task-emotion
settings.firstlevel.emotion.run1.sess_params.timing_units = 'secs';
settings.firstlevel.emotion.run1.sess_params.timing_RT = 2;
settings.firstlevel.emotion.run1.sess_params.cond_names = {'Faces', 'Shapes'};
settings.firstlevel.emotion.run2.sess_params.timing_units = 'secs';
settings.firstlevel.emotion.run2.sess_params.timing_RT = 2;
settings.firstlevel.emotion.run2.sess_params.cond_names = {'MentalEmotion'};

% Settings for plotting task conditions
onset = [11; 31; 51; 71; 91; 111; 131; 151; 171; 191];
duration = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10];
settings.firstlevel.motor.run1.plot_params.cond_onset = onset;
settings.firstlevel.motor.run1.plot_params.cond_duration = duration;
settings.firstlevel.motor.run2.plot_params.cond_onset = onset;
settings.firstlevel.motor.run2.plot_params.cond_duration = duration;
settings.firstlevel.emotion.run2.plot_params.cond_onset = onset;
settings.firstlevel.emotion.run2.plot_params.cond_duration = duration;
onset = [12; 32; 52; 72; 92; 112; 132; 152; 172; 192];
duration = [9; 9; 9; 9; 9; 9; 9; 9; 9; 9];
settings.firstlevel.emotion.run1.plot_params.cond_onset = onset;
settings.firstlevel.emotion.run1.plot_params.cond_duration = duration;

% Settings for first level analysis: glm regressors to include
settings.firstlevel.glm_regressors.trans_rot = true;
settings.firstlevel.glm_regressors.trans_rot_derivative1 = true;
settings.firstlevel.glm_regressors.trans_rot_power2 = false;
settings.firstlevel.glm_regressors.trans_rot_derivative1_power2 = false;
settings.firstlevel.glm_regressors.framewise_displacement_censor02 = false;
settings.firstlevel.glm_regressors.framewise_displacement_censor05 = false;
settings.firstlevel.glm_regressors.dvars_censor = false; % not yet implemented
settings.firstlevel.glm_regressors.std_dvars_censor = false; % not yet implemented
settings.firstlevel.glm_regressors.grey_matter = false;
settings.firstlevel.glm_regressors.white_matter = false;
settings.firstlevel.glm_regressors.csf = true;
settings.firstlevel.glm_regressors.global_signal = false;
% Order of included retroicor regressors; if 0 ==> exclude
settings.firstlevel.glm_regressors.retroicor_c = 2; % cardiac, max 6
settings.firstlevel.glm_regressors.retroicor_r = 2; % respiratory, max 8
settings.firstlevel.glm_regressors.retroicor_cxr = 0; % interaction, max 4
settings.firstlevel.glm_regressors.hrv = false;
settings.firstlevel.glm_regressors.rvt = false;


% Settings for first level analysis: task-motor
settings.firstlevel.motor.run1.contrast_params.consess{1}.tcon.name = 'FingerTapping';
settings.firstlevel.motor.run1.contrast_params.consess{1}.tcon.weights = [1];
settings.firstlevel.motor.run1.contrast_params.consess{1}.tcon.sessrep = 'none';
settings.firstlevel.motor.run2.contrast_params.consess{1}.tcon.name = 'MentalFingerTapping';
settings.firstlevel.motor.run2.contrast_params.consess{1}.tcon.weights = [1];
settings.firstlevel.motor.run2.contrast_params.consess{1}.tcon.sessrep = 'none';

% Settings for first level analysis: task-emotion
settings.firstlevel.emotion.run1.contrast_params.consess{1}.tcon.name = 'Faces';
settings.firstlevel.emotion.run1.contrast_params.consess{1}.tcon.weights = [1 0];
settings.firstlevel.emotion.run1.contrast_params.consess{1}.tcon.sessrep = 'none';
settings.firstlevel.emotion.run1.contrast_params.consess{2}.tcon.name = 'Shapes';
settings.firstlevel.emotion.run1.contrast_params.consess{2}.tcon.weights = [0 1];
settings.firstlevel.emotion.run1.contrast_params.consess{2}.tcon.sessrep = 'none';
settings.firstlevel.emotion.run1.contrast_params.consess{3}.tcon.name = 'Faces>Shapes';
settings.firstlevel.emotion.run1.contrast_params.consess{3}.tcon.weights = [1 -1];
settings.firstlevel.emotion.run1.contrast_params.consess{3}.tcon.sessrep = 'none';
settings.firstlevel.emotion.run2.contrast_params.consess{1}.tcon.name = 'MentalEmotion';
settings.firstlevel.emotion.run2.contrast_params.consess{1}.tcon.weights = [1];
settings.firstlevel.emotion.run2.contrast_params.consess{1}.tcon.sessrep = 'none';

%matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Patients > Control';
%matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [-1 1];
%matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';