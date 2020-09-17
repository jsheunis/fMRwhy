% A custom workflow that runs 1st level analysis for all runs of all tasks of specified subjects

% Code steps:
% 1.


%--------------------------------------------------------------------------


% -------
% STEP 0.1 -- Load defaults, filenames and parameters
% -------

% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Loop through subjects, sessions, tasks, runs, etc
%subs = {'001', '002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015'};
subs = {'001', '002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
subs = {'019'};
ses = '';
%tasks = {'motor', 'emotion'};
%runs = {'1', '2'};
%echoes = {'2', 'combinedMEtsnr', 'combinedMEt2star', 'combinedMEte'};


for s = 1:numel(subs)

    sub = subs{s};
    % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

    fs_sub_dir = fullfile(options.deriv_dir, 'freesurfer', ['sub-' sub]);

    % to get template file
    ses = '';
    task = 'rest';
    run = '1';
    echo = '2';
    options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

    source_fn = options.anatomical_fn;
    reference_fn = options.template_fn;
    fn1 = fullfile(fs_sub_dir, 'BA_4a_lh.nii');
    fn2 = fullfile(fs_sub_dir, 'BA_4p_lh.nii');
    other_fn = {fn1, fn2};
    interp = 0; % nearest neighbour
    saveAs_fn = '';

    fmrwhy_batch_coregEstResl(source_fn, reference_fn, other_fn, interp, saveAs_fn)

end