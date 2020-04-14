function output = fmrwhy_tapas_PhysIO(options)

%options.save_dir = '';
%options.vendor = 'BIDS';
%options.cardiac_fn = '';
%options.respiration_fn = '';
%options.sampling_interval = 0.002; % 500 Hz ==> Philips wired acquisition
%options.align_scan = 'last';
%options.Nslices = 34;
%options.TR = 2; % in seconds
%options.Ndummies = 5; % include, even if these are not included in the fMRI timeseries data exported from the scanner
%options.Nscans = 210;
%options.onset_slice = 1;
%options.cardiac_modality = 'PPU';
%options.output_multiple_regressors_fn = ''; % text file name


physio = tapas_physio_new('empty');


physio.save_dir = {options.save_dir};
physio.log_files.vendor = options.vendor;
physio.log_files.cardiac = {options.cardiac_fn};
physio.log_files.respiration = {options.respiration_fn};
physio.log_files.scan_timing = {''};
physio.log_files.sampling_interval = options.sampling_interval;
physio.log_files.relative_start_acquisition = 0;
physio.log_files.align_scan = options.align_scan;
physio.scan_timing.sqpar.Nslices = options.Nslices;
physio.scan_timing.sqpar.NslicesPerBeat = [];
physio.scan_timing.sqpar.TR = options.TR;
physio.scan_timing.sqpar.Ndummies = options.Ndummies;
physio.scan_timing.sqpar.Nscans = options.Nscans;
physio.scan_timing.sqpar.onset_slice = options.onset_slice;
physio.scan_timing.sqpar.time_slice_to_slice = [];
physio.scan_timing.sqpar.Nprep = [];
physio.scan_timing.sync.nominal = [];
% physio.scan_timing.sync.gradient_log.grad_direction = 'y';
% physio.scan_timing.sync.gradient_log.zero = 700;
% physio.scan_timing.sync.gradient_log.slice = 1800;
% physio.scan_timing.sync.gradient_log.vol = [];
% physio.scan_timing.sync.gradient_log.vol_spacing = [];
physio.preproc.cardiac.modality = 'PPU';
physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
physio.preproc.cardiac.initial_cpulse_select.auto_matched.max_heart_rate_bpm = 90;
physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
physio.model.output_multiple_regressors = options.output_multiple_regressors_fn;
physio.model.output_physio = 'physio.mat';
physio.model.orthogonalise = 'none';
physio.model.censor_unreliable_recording_intervals = false;
physio.model.retroicor.yes.order.c = 3;
physio.model.retroicor.yes.order.r = 4;
physio.model.retroicor.yes.order.cr = 1;
physio.model.rvt.yes.delays = 0;
physio.model.hrv.yes.delays = 0;
physio.model.noise_rois.no = struct([]);
physio.model.movement.no = struct([]);
physio.model.other.no = struct([]);
physio.verbose.level = options.level;
physio.verbose.fig_output_file = options.fig_output_file;
physio.verbose.use_tabs = false;

[physio, R, ons_secs] = tapas_physio_main_create_regressors(physio)

output = struct;
output.physio_mat_fn = fullfile(options.save_dir, 'physio.mat');
output.physio_mat = physio;
%output.physio_mat = load(output.physio_mat_fn);
output.multiple_regressors_fn = fullfile(options.save_dir, options.output_multiple_regressors_fn);
output.physio_regr = R;
%output.physio_regr = load(output.multiple_regressors_fn);


% Regressors in (output_multiple_regressors_fn)
%RETROICOR cardiac regressors [2 x nOrderCardiac] ==> 2 x 3 = 6
%RETROICOR respiratory regressors [2 x nOrderRespiratory] ==> 2 x 4 = 8
%RETROICOR cardXResp interaction regressors [4 x nOrderCardiacXRespiratory] ==> 4 x 1 = 4
%HRV [nDelaysHRV]
%RVT [nDelaysRVT]
%Noise ROIs (PCA signatures and mean of each region) [nNoiseROIs x (nComponents+1)]
%Other (included other text file) [nColumnsOtherFile]
%Motion [6 or 12 or 24, depending on motion model]