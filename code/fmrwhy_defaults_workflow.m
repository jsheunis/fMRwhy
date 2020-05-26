function defaults = fmrwhy_defaults_workflow(bids_dir)

% TODO: this file needs rework. Many of the parameters are specific to a single workflow.
% and many others are standard for fMRwhy. These need to be separated.

% SPM directory
defaults.spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';

% Derivatives directory setup
defaults.deriv_dir = fullfile(bids_dir, 'derivatives');
defaults.preproc_dir = fullfile(defaults.deriv_dir, 'fmrwhy-preproc');
defaults.qc_dir = fullfile(defaults.deriv_dir, 'fmrwhy-qc');
defaults.stats_dir = fullfile(defaults.deriv_dir, 'fmrwhy-stats');

% Set template for functional realignment purposes
defaults.template_task = 'rest';
defaults.template_run = '1';
defaults.template_echo = '2';

% Sequence parameters
defaults.TR = 2;
defaults.N_slices = 34;
defaults.Ndummies = 5
defaults.Nscans = 210

% Dataset parameters
defaults.Nsessions = 2;
defaults.Nruns = 2;
defaults.tasks = {'rest', 'motor', 'emotion'};
defaults.Ntasks = numel(defaults.tasks);

% Settings for structFunc processing
defaults.has_rois = 1;
defaults.roi_orig_dir = '/Volumes/Stephan_WD/NEUFEPME_data_templates';
defaults.roi_orig_fns = {'',
                   {'Left_Motor_4a_4p.nii', 'Right_Motor_4a_4p.nii'},
                   {'Left_Amygdala_allregions.nii', 'Right_Amygdala_allregions.nii'}} % cell array of cell arrays, in same order as tasks
defaults.roi_names = {'',
                   {'Left Motor', 'Right Motor'},
                   {'Left Amygdala', 'Right Amygdala'}} % cell array of cell arrays, in same order as tasks
defaults.roi_desc = {'',
                   {'leftMotor', 'rightMotor'},
                   {'leftAmygdala', 'rightAmygdala'}} % cell array of cell arrays, in same order as tasks

%defaults.roi_fns = ROIs in subject space (not resliced)
%defaults.rroi_fns = resliced ROIs in subject space


% Settings for basicFunc processing
defaults.fwhm = 7;

% Settings for generateMultRegr routine
defaults.confounds.include_fd = 1;
defaults.confounds.include_tissue = 1;
defaults.confounds.include_physio = 0;

% generateMultRegr: framewise displacement
defaults.r = 50; % mm
defaults.FD_threshold = 0.5; % mm

% generateMultRegr: PhysIO
defaults.physio.options.save_dir = '';
defaults.physio.options.cardiac_fn = '';
defaults.physio.options.respiration_fn = '';
defaults.physio.options.vendor = 'BIDS';
defaults.physio.options.sampling_interval = 0.002; % 500 Hz ==> Philips wired acquisition
defaults.physio.options.align_scan = 'last';
defaults.physio.options.Nslices = defaults.N_slices;
defaults.physio.options.TR = defaults.TR; % in seconds
defaults.physio.options.Ndummies = defaults.Ndummies; % include, even if these are not included in the fMRI timeseries data exported from the scanner
defaults.physio.options.Nscans = defaults.Nscans;
defaults.physio.options.onset_slice = 1;
defaults.physio.options.cardiac_modality = 'PPU';
defaults.output_multiple_regressors_fn = 'PhysIO_multiple_regressors.txt'; % text file name

% Settings for QC
defaults.theplot.intensity_scale = [-6 6];



%% Settings for first level analysis: task-motor
%defaults.tasks.motor.sess_params = struct;
%defaults.tasks.motor.sess_params.timing_units = 'scans';
%defaults.tasks.motor.sess_params.timing_RT = 2;
%defaults.tasks.motor.sess_params.cond_name = 'Finger_tapping_rhs';
%defaults.tasks.motor.sess_params.cond_onset = [11; 31; 51; 71; 91; 111; 131; 151; 171; 191];
%defaults.tasks.motor.sess_params.cond_duration = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10];
%
%% Settings for first level analysis: task-emotion
%defaults.tasks.emotion.sess_params = struct;
%defaults.tasks.emotion.sess_params.timing_units = 'scans';
%defaults.tasks.emotion.sess_params.timing_RT = 2;
%defaults.tasks.emotion.sess_params.cond_name = 'Match_shapes_faces';
%defaults.tasks.emotion.sess_params.cond_onset = [11; 31; 51; 71; 91; 111; 131; 151; 171; 191];
%defaults.tasks.emotion.sess_params.cond_duration = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10];
