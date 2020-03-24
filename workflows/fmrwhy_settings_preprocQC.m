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
settings.template_run = '1';
settings.template_echo = '2';

% Sequence parameters
settings.TR = 2;
settings.N_slices = 34;
settings.Ndummies = 5
settings.Nscans = 210

% Dataset parameters
settings.Nsessions = 2;
settings.Nruns = 2;
settings.tasks = {'rest', 'motor', 'emotion'};
settings.Ntasks = numel(settings.tasks);

% Settings for structFunc processing

% Settings for anatLocaliser processing
settings.map_rois = 1;
settings.roi_orig_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_templates';
settings.roi = struct;
settings.roi.motor.orig_fn = {fullfile(settings.roi_orig_dir, 'Left_Motor_4a_4p.nii'),
                                fullfile(settings.roi_orig_dir, 'Right_Motor_4a_4p.nii')}; % Raw ROI filenames

settings.roi.motor.name = {'Left Motor', 'Right Motor'}; % For plots and strings
settings.roi.motor.desc = {'leftMotor', 'rightMotor'}; % For BIDS file naming (after normalisation  to functional space)

settings.roi.emotion.orig_fn = {fullfile(settings.roi_orig_dir, 'Left_Motor_4a_4p.nii'),
                                fullfile(settings.roi_orig_dir, 'Right_Motor_4a_4p.nii')}; % Raw ROI filenames

settings.roi.emotion.name = {'Left Amygdala', 'Right Amygdala'}; % For plots and strings
settings.roi.emotion.desc = {'leftAmygdala', 'rightAmygdala'}; % For BIDS file naming (after normalisation  to functional space)

%settings.roi.(task).roi_fn = ROIs in subject space (not resliced)
%settings.roi.(task).rroi_fn = resliced ROIs in subject space


% Settings for basicFunc processing
settings.fwhm = 7;

% Settings for generateMultRegr routine
settings.confounds.include_fd = 1;
settings.confounds.include_tissue = 1;
settings.confounds.include_physio = 0;

% generateMultRegr: framewise displacement
settings.r = 50; % mm
settings.FD_threshold = 0.5; % mm

% generateMultRegr: PhysIO
settings.physio.options.save_dir = '';
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

% Settings for QC
settings.theplot.intensity_scale = [-6 6];



%% Settings for first level analysis: task-motor
%settings.tasks.motor.sess_params = struct;
%settings.tasks.motor.sess_params.timing_units = 'scans';
%settings.tasks.motor.sess_params.timing_RT = 2;
%settings.tasks.motor.sess_params.cond_name = 'Finger_tapping_rhs';
%settings.tasks.motor.sess_params.cond_onset = [11; 31; 51; 71; 91; 111; 131; 151; 171; 191];
%settings.tasks.motor.sess_params.cond_duration = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10];
%
%% Settings for first level analysis: task-emotion
%settings.tasks.emotion.sess_params = struct;
%settings.tasks.emotion.sess_params.timing_units = 'scans';
%settings.tasks.emotion.sess_params.timing_RT = 2;
%settings.tasks.emotion.sess_params.cond_name = 'Match_shapes_faces';
%settings.tasks.emotion.sess_params.cond_onset = [11; 31; 51; 71; 91; 111; 131; 151; 171; 191];
%settings.tasks.emotion.sess_params.cond_duration = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10];
