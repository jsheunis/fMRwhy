% A custom workflow that does ...


%--------------------------------------------------------------------------

% Load/create required parameters
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Loop through subjects, sessions, tasks, runs, etc
sub = '001';
ses = '';
task = 'motor';
%task = 'emotion';
run = '1';
echo = '2';

options = struct;

fmrwhy_workflow_1stlevel(bids_dir, sub, ses, task, run, echo, options)