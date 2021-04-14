% Load/create required parameters
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

%% Setup fmrwhy BIDS-derivatuve directories on workflow level
% options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);
%
%% Grab parameters from workflow settings file
% options = fmrwhy_settings_preprocQC(bids_dir, options);

% Loop through subjects, sessions, tasks, runs, etc
sub = '001';

%% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
% options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);
%
%% Update workflow params with subject anatomical derivative filenames
% options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% Loop through sessions, tasks, runs, etc
ses = '';
task = 'rest';
run = '1';
echo = '2';

%% Update workflow params with subject functional derivative filenames
% options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

fmrwhy_qc_createMaskMontages(bids_dir, sub, 1);
