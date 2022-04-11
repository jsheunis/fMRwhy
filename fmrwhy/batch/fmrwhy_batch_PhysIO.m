function output = fmrwhy_batch_PhysIO(options)

    % options.save_dir = '';
    % options.vendor = 'BIDS';
    % options.cardiac_fn = '';
    % options.respiration_fn = '';
    % options.sampling_interval = 0.002; % 500 Hz ==> Philips wired acquisition
    % options.align_scan = 'last';
    % options.Nslices = 34;
    % options.TR = 2; % in seconds
    % options.Ndummies = 5; % include, even if these are not included in the fMRI timeseries data exported from the scanner
    % options.Nscans = 210;
    % options.onset_slice = 1;
    % options.cardiac_modality = 'PPU';
    % options.output_multiple_regressors_fn = ''; % text file name

    spm('defaults', 'fmri');
    spm_jobman('initcfg');

    phys = struct;
    phys.matlabbatch{1}.spm.tools.physio.save_dir = {options.save_dir};
    phys.matlabbatch{1}.spm.tools.physio.log_files.vendor = options.vendor;
    phys.matlabbatch{1}.spm.tools.physio.log_files.cardiac = {options.cardiac_fn};
    phys.matlabbatch{1}.spm.tools.physio.log_files.respiration = {options.respiration_fn};
    phys.matlabbatch{1}.spm.tools.physio.log_files.scan_timing = {''};
    phys.matlabbatch{1}.spm.tools.physio.log_files.sampling_interval = options.sampling_interval;
    phys.matlabbatch{1}.spm.tools.physio.log_files.relative_start_acquisition = 0;
    phys.matlabbatch{1}.spm.tools.physio.log_files.align_scan = options.align_scan;
    phys.matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nslices = options.Nslices;
    phys.matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
    phys.matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.TR = options.TR;
    phys.matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Ndummies = options.Ndummies;
    phys.matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nscans = options.Nscans;
    phys.matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.onset_slice = options.onset_slice;
    phys.matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
    phys.matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
    phys.matlabbatch{1}.spm.tools.physio.scan_timing.sync.nominal = [];
    % phys.matlabbatch{1}.spm.tools.physio.scan_timing.sync.gradient_log.grad_direction = 'y';
    % phys.matlabbatch{1}.spm.tools.physio.scan_timing.sync.gradient_log.zero = 700;
    % phys.matlabbatch{1}.spm.tools.physio.scan_timing.sync.gradient_log.slice = 1800;
    % phys.matlabbatch{1}.spm.tools.physio.scan_timing.sync.gradient_log.vol = [];
    % phys.matlabbatch{1}.spm.tools.physio.scan_timing.sync.gradient_log.vol_spacing = [];
    
    phys.matlabbatch{1}.spm.tools.physio.preproc.cardiac.modality = options.cardiac_modality;
    phys.matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
    phys.matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
    phys.matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.max_heart_rate_bpm = 90;
    phys.matlabbatch{1}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
    phys.matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = options.output_multiple_regressors_fn;
    phys.matlabbatch{1}.spm.tools.physio.model.output_physio = 'physio.mat';
    phys.matlabbatch{1}.spm.tools.physio.model.orthogonalise = 'none';
    phys.matlabbatch{1}.spm.tools.physio.model.censor_unreliable_recording_intervals = false;
    phys.matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.c = 3;
    phys.matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.r = 4;
    phys.matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.cr = 1;
    phys.matlabbatch{1}.spm.tools.physio.model.rvt.yes.delays = 0;
    phys.matlabbatch{1}.spm.tools.physio.model.hrv.yes.delays = 0;
    phys.matlabbatch{1}.spm.tools.physio.model.noise_rois.no = struct([]);
    phys.matlabbatch{1}.spm.tools.physio.model.movement.no = struct([]);
    phys.matlabbatch{1}.spm.tools.physio.model.other.no = struct([]);
    phys.matlabbatch{1}.spm.tools.physio.verbose.level = options.level;
    phys.matlabbatch{1}.spm.tools.physio.verbose.fig_output_file = options.fig_output_file;
    phys.matlabbatch{1}.spm.tools.physio.verbose.use_tabs = false;

    spm_jobman('run', phys.matlabbatch);
    output = struct;
    output.physio_mat_fn = fullfile(options.save_dir, 'physio.mat');
    output.physio_mat = load(output.physio_mat_fn);
    output.multiple_regressors_fn = fullfile(options.save_dir, options.output_multiple_regressors_fn);
    output.physio_regr = load(output.multiple_regressors_fn);

    % Regressors in (output_multiple_regressors_fn)
    % RETROICOR cardiac regressors [2 x nOrderCardiac] ==> 2 x 3 = 6
    % RETROICOR respiratory regressors [2 x nOrderRespiratory] ==> 2 x 4 = 8
    % RETROICOR cardXResp interaction regressors [4 x nOrderCardiacXRespiratory] ==> 4 x 1 = 4
    % HRV [nDelaysHRV]
    % RVT [nDelaysRVT]
    % Noise ROIs (PCA signatures and mean of each region) [nNoiseROIs x (nComponents+1)]
    % Other (included other text file) [nColumnsOtherFile]
    % Motion [6 or 12 or 24, depending on motion model]
