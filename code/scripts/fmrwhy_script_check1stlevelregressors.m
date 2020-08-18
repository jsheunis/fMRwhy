% A custom workflow that runs 1st level analysis for all runs of all tasks of specified subjects

% Code steps:
% 1.


%--------------------------------------------------------------------------


% -------
% STEP 0.1 -- Load defaults, filenames and parameters
% -------

% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Loop through subjects, sessions, tasks, runs, etc
subs = {'001', '002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015'};
%subs = {'001', '002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
%sub = '002';
ses = '';
tasks = {'motor', 'emotion'};
runs = {'1', '2'};
echoes = {'2', 'combinedMEtsnr', 'combinedMEt2star', 'combinedMEte'};


for s = 1:numel(subs)
    sub = subs{s};
    % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

    % -------
    % PER TASK and RUN
    % -------
    % Loop through sessions, tasks, runs, echoes.

    toTransform_fns = {};
    saveAs_transform_fns = {};
    for t = 1:numel(tasks)
        task = tasks{t};
        for r = 1:numel(runs)
            run = runs{r};
            for e = 1:numel(echoes)
                echo = echoes{e};

                options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, options.template_echo, options);
                options.sub_dir_stats = fullfile(options.stats_dir, ['sub-' sub]);
                run_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-' echo]);

                regr_fn = fullfile(run_dir_stats,['sub-' sub '_task-' task '_run-' run '_echo-2_desc-GLM_regressors.txt']);

                regressors = load(regr_fn);
                [x, y] = size(regressors);

                if y~= 17
                    txt = ['sub-' sub '_task-' task '_run-' run '_echo-' echo ': ' num2str(y)];
                    disp(txt);
                end

            end
        end
    end
end