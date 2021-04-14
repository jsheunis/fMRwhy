function fmrwhy_qc_customWF(bids_dir, sub, ses, task, run, template)

    % Software dependences:
    % - Matlab vX
    % - SPM12 r7771
    % - (TAPAS PhysIO Toolbox vX)

    % Data required in order for this function to run as intended:
    % - T1w (raw data)
    % - Single run of fMRI timeseries (raw data)
    % - Head motion/movement parameters derived from unprocessed data (i.e. from realignment of raw fMRI timeseries)
    % - T1w coregistered to template functional volume (from fmrwhy_preproc_structFunc.m))
    % - Segmentations of GM, WM, CSF in template functional volume space (from fmrwhy_preproc_structFunc.m))

    % Pseudo-code:
    % 1. Get BIDS directory of run
    % 2. Create relevant derivative directories (qc and preproc) if it doesn't exist yet
    % 3. Run minimal preproc on run for Qc:
    %   - Head motion/movement parameters derived from unprocessed data (i.e. from realignment of raw fMRI timeseries)
    %   - Framewise displacement
    %   - Statistical measures / images
    %   -
    %   -
    %   -
    % 4.

    % Input params
    bids_dir = '/Volumes/Stephan_WD/NEUFEPME_data_BIDS';
    sub = '001';
    ses = '';
    task = '';
    run = 1;
    echo = 2;

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
    if ~exist(sub_dir_preproc, 'dir')
        mkdir(sub_dir_preproc);
        sub_dir_BIDS = fullfile(bids_dir, ['sub-' sub]);
        copyfile(sub_dir_BIDS, sub_dir_preproc);
    end
    if ~exist(sub_dir_qc, 'dir')
        mkdir(sub_dir);
    end

    % spm_BIDS(BIDS,'data','sub','005','task','rest','echo','2','type','bold')
