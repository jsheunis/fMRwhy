% Load/create required parameters
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Setup fmrwhy bids directories on workflow level
fmrwhy_defaults_setupDerivDirs(bids_dir);

% Grab default workflow params
wf_params = fmrwhy_defaults_workflow(bids_dir);

% Loop through subjects, sessions, tasks, runs, etc
sub = '001';

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
fmrwhy_defaults_setupSubDirs(bids_dir, sub);

% Update workflow params with subject anatomical derivative filenames
wf_params = fmrwhy_defaults_subAnat(bids_dir, sub, wf_params);

% Loop through sessions, tasks, runs, etc
ses = '';
task = 'motor';
run = '1';
echo = '2';

% Update workflow params with subject functional derivative filenames
wf_params = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, wf_params);

fmrwhy_qc_createMaskMontages(bids_dir, sub, 0)