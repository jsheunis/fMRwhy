
options = fmrwhy_defaults;
% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';
% Loop through subjects, sessions, tasks, runs, etc
sub = '001';
% Loop through sessions, tasks, runs, etc
ses = '';


% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% TODO: have to loop through runs for this one if we want to generate report for all runs for sub:
% -------
% PER TASK and RUN
% -------

% Loop through sessions, tasks, runs, echoes.
tasks = {'rest', 'motor', 'emotion'};
runs = {'1', '2'};
% Note: Choose arbitrary echo for now, since this is not needed for current qc report
echo = options.template_echo;

for t = 1:numel(tasks)
    task = tasks{t};
    for r = 1:numel(runs)
        run = runs{r};
        % Update workflow params with subject functional derivative filenames
        options = fmrwhy_defaults_subFunc(bids_dir, sub, '', task, run, echo, options);
        fmrwhy_qc_generateSubRunReport(bids_dir, sub, task, run, options)
    end
end

