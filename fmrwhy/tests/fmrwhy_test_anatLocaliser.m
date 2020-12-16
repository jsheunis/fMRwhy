% A custom workflow that does ...


%--------------------------------------------------------------------------

% Load/create required parameters
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Loop through subjects, sessions, tasks, runs, etc
sub = '001';

options = struct;

fmrwhy_preproc_anatLocaliser(bids_dir, sub, options)