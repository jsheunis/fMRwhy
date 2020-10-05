function options = fmrwhy_settings_validate(options)
% fmrwhy_settings_validate: Validate and derive settings for the fmrwhy_workflow_qc pipeline

% BIDS structure values
options.bids_dataset = bids.layout(options.bids_dir);
options.subjects = bids.query(options.bids_dataset,'subjects');
options.sessions = bids.query(options.bids_dataset,'sessions');
options.runs = bids.query(options.bids_dataset,'runs');
options.tasks = bids.query(options.bids_dataset,'tasks');
options.types = bids.query(options.bids_dataset,'types');
options.modalities = bids.query(options.bids_dataset,'modalities');

% set list of subjects for which to generate QC output
% If options.subjects_output is set to 'all' (default), read the list from the bids.query output, otherwise leave as is 
if strcmp(options.subjects_output, 'all')
    options.subjects_output = options.subjects;
end

% Set template for functional realignment purposes (if not needed, set to [])
options.template_sub = options.subjects{1};

% Dataset parameters
options.Nsessions = numel(options.sessions);
options.Ntasks = numel(options.tasks);
options.Nruns = numel(options.runs);

% Derive template flags from BIDS structure
options.has_sessions = ~isempty(options.sessions);
options.has_runs = ~isempty(options.runs);
if options.has_sessions
    if options.has_runs
        filenames = bids.query(options.bids_dataset, 'data', 'sub', options.template_sub, 'task', options.template_task, 'sess', options.template_session, 'run', options.template_run, 'type', 'bold');
    else
        filenames = bids.query(options.bids_dataset, 'data', 'sub', options.template_sub, 'task', options.template_task, 'sess', options.template_session, 'type', 'bold');
    end
else
    if options.has_runs
        filenames = bids.query(options.bids_dataset, 'data', 'sub', options.template_sub, 'task', options.template_task, 'run', options.template_run, 'type', 'bold');
    else
        filenames = bids.query(options.bids_dataset, 'data', 'sub', options.template_sub, 'task', options.template_task, 'type', 'bold');
    end
end
options.N_echoes = numel(filenames);
options.is_multiecho = false;
if options.N_echoes > 1
    options.is_multiecho = true;
end

% Sequence parameters
options.Ne = numel(options.TE);

if options.is_multiecho
    if options.N_echoes ~= options.Ne
        disp('ERROR: number of echoes derived from BIDS dataset (using bids-matlab) do not match the number of echo times specified in settings file. FIX!')
    end
end