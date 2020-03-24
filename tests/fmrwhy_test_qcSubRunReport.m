
options = fmrwhy_defaults;
% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';
% Loop through subjects, sessions, tasks, runs, etc
sub = '001';
% Loop through sessions, tasks, runs, etc
ses = '';
task = 'rest';
run = '1';
echo = '2';

fmrwhy_qc_generateSubRunReport(bids_dir, sub, task, run, options)