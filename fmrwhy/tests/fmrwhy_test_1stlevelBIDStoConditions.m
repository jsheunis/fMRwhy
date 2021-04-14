% A custom workflow that does ...

% --------------------------------------------------------------------------

% Load/create required parameters
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Loop through subjects, sessions, tasks, runs, etc
sub = '001';
ses = '';
task = 'emotion';
run = '1';
echo = '2';

events_fn = fullfile(bids_dir, 'derivatives', 'fmrwhy-preproc', ['sub-' sub], 'func',  ['sub-' sub '_task-' task '_run-' run '_events.tsv']);
cond_names = {'Shapes', 'Faces'};

[cond, trials] = fmrwhy_util_1stlevelBIDStoConditions(events_fn, cond_names);
