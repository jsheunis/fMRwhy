function output = fmrwhy_qc_generatePhysRegr(bids_dir, sub, ses, task, run)

    % BIDS structure values
    BIDS = spm_BIDS(bids_dir);
    subjects = spm_BIDS(BIDS, 'subjects');
    sessions = spm_BIDS(BIDS, 'sessions');
    runs = spm_BIDS(BIDS, 'runs');
    tasks = spm_BIDS(BIDS, 'tasks');
    types = spm_BIDS(BIDS, 'types');
    modalities = spm_BIDS(BIDS, 'modalities');

    % Directory and content setup
    deriv_dir = fullfile(bids_dir, 'derivatives');
    preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
    qc_dir = fullfile(deriv_dir, 'fmrwhy-qc');
    sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);
    sub_dir_qc = fullfile(qc_dir, ['sub-' sub]);
    func_dir_qc = fullfile(sub_dir_qc, 'func');
    if ~exist(func_dir_qc, 'dir')
        mkdir(func_dir_qc);
    end

    % Load physio log-file, e.g. sub-001_task-rest_run-1_physio.tsv.gz
    log_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' num2str(run) '_physio.tsv.gz']);

    % Create struct with options for PhysIO, ideally from
    options.save_dir = func_dir_qc;
    options.vendor = 'BIDS';
    options.cardiac_fn = log_fn;
    options.respiration_fn = log_fn;
    options.sampling_interval = 0.002; % 500 Hz ==> Philips wired acquisition
    options.align_scan = 'last';
    options.Nslices = 34;
    options.TR = 2; % in seconds
    options.Ndummies = 5; % include, even if these are not included in the fMRI timeseries data exported from the scanner
    options.Nscans = 210;
    options.onset_slice = 1;
    options.cardiac_modality = 'PPU';
    options.output_multiple_regressors_fn = 'PhysIO_multiple_regressors.txt'; % text file name

    % Run PhysIO
    phys_data = fmrwhy_batch_PhysIO(options);
    col_names = {'retroicor_c1', 'retroicor_c2', 'retroicor_c3', 'retroicor_c4', 'retroicor_c5', 'retroicor_c6', 'retroicor_r1', 'retroicor_r2', 'retroicor_r3', 'retroicor_r4', 'retroicor_r5', 'retroicor_r6', 'retroicor_r7', 'retroicor_r8', 'retroicor_cxr1', 'retroicor_cxr2', 'retroicor_cxr3', 'retroicor_cxr4', 'hrv', 'rvt'};
    output.multiple_regressors_fn = fmrwhy_util_saveAsTSV(phys_data.multiple_regressors_fn, col_names);
