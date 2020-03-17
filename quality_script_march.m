bids_dir = '/Volumes/Stephan_WD/NEUFEPME_data_BIDS';

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
task = 'rest';
run = '1';
echo = '2';

% Update workflow params with subject functional derivative filenames
wf_params = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, wf_params);

%%
% -------
% STEP 2: Calculate and generate statistical measures / images (tsnr, variance, std, psc, DVARS)
% -------
stats = fmrwhy_qc_createStatsOutput(bids_dir, sub, ses, task, run, echo, wf_params);
%%
% -------
% STEP 3: Plot The Plot
% -------
fmrwhy_qc_createThePlot(bids_dir, sub, ses, task, run, echo, wf_params);