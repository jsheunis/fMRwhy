% A custom workflow that does structFunc and basicFunc preprocessing and QC for a single subject in the NEUFEP study

% Code steps:
% 1. Define template/default variables, directories and filenames
% 2. Specify


%--------------------------------------------------------------------------


% -------
% STEP 0.1 -- Load defaults, filenames and parameters
% -------

% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Loop through subjects, sessions, tasks, runs, etc
sub = '001';

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% Grab functional template
template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
options.template_fn = template_fn;

% -------
% 1ST LEVEL ANALYSIS PER TASK AND RUN
% -------

% Loop through sessions, tasks, runs, etc
ses = '';
tasks = {'rest', 'motor', 'emotion'};


for t = 1:numel(tasks)

    task = tasks{t};
    run = '1';
    echo = '2';

    % Update workflow params with subject functional derivative filenames
    options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

    % --------------------------------------------------------------------------
    % STEP 1 -- Structural-functional preprocessing: fmrwhy_preproc_structFunc.m
    % --------------------------------------------------------------------------





    %
    %% Step 3: anatomical-localizer-preproc:     - rtme_preproc_anatLocaliser.m
    %fmrwhy_preproc_anatLocaliser(sub, options);

    % Step 3: functional-localizer-preproc:     - rtme_preproc_funcLocaliser.m
    %                                           - rtme_preproc_generateRegressors.m
    %                                           - rtme_preproc_generateRetroicor.m
    %                                           - rtme_preproc_generateFDregr.m
    %                                           - rtme_preproc_generateTissueSignals.m

    %
    %for t = 1:numel(defaults.tasks)
    %%    disp(['Performing 3D volume realignment for: ' sub '_task-' tasks(t) '_run-' defaults.template_run])
    %    rtme_preproc_funcLocaliser(sub, task, defaults.template_run, options.template_echo, defaults)
    %end


    % Step 4: calculate-prior-measures-preproc - rtme_preproc_estimateParams.m
    %rtme_preproc_estimateParams(sub, defaults);
end



