function fmrwhy_qc_sub(bids_dir, sub, ses, task, run, echo)
    % A pipeline to run QC processes and generate metrics/visualisations for a single subject,
    % assuming that there are

    % Software dependences:
    % - Matlab vX
    % - SPM12 r7771
    % - (TAPAS PhysIO Toolbox vX)

    % Data required in order for this function to run as intended:
    % - T1w (raw data)
    % - fMRI timeseries (raw data)
    % - T1w coregistered to template functional volume (from fmrwhy_preproc_structFunc.m)
    % - Segmentations of GM, WM, CSF in template functional volume space (from fmrwhy_preproc_structFunc.m)
    % - Head motion/movement parameters derived from unprocessed data (i.e. from realignment of raw fMRI timeseries)

    % Code steps:
    % 1. QC steps for structural-functional preproc:
    %   - registration
    %   - what else?

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

    % spm_BIDS(BIDS,'data','sub','005','task','rest','echo','2','type','bold')

    %%

    % Directory and content setup
    deriv_dir = fullfile(bids_dir, 'derivatives');
    preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
    qc_dir = fullfile(deriv_dir, 'fmrwhy-qc');
    sub_dir_qc = fullfile(qc_dir, ['sub-' sub]);

    if ~exist(sub_dir_qc, 'dir')
        mkdir(sub_dir_qc);
    end

    % 1 FIRST DO STRUCTFUNC QC
    %   - registration: mask contour overlaid on template EPI, or mean EPI
    %   - others? look at QPAC, MRIQC, etc

    % Plot masks on montage of template EPI to show struct-func registration
    montage_data = fmrwhy_qc_createMaskMontages(bids_dir, sub, 1);
    % Save images to relevant directories

    % 2 THEN DO fmrwhy_qc_run FOR EACH RUN IN SUBJECT FOLDER
    % for i = 1:numel(runs)
    ses = '';
    task = '';
    run = 1;
    echo = 2;
    fmrwhy_qc_run(bids_dir, sub, ses, task, run, template);
    % end
