% A custom workflow that runs 2nd level analysis

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
subs = {'001', '002', '003', '004', '005'};
%subs = {'001', '002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
%sub = '002';
ses = '';
tasks = {'motor'};
%tasks = {'motor', 'emotion'};
runs = {'1'};
%runs = {'1', '2'};
echoes = {'2'};
%echoes = {'2', 'combinedMEtsnr', 'combinedMEt2star', 'combinedMEte'};

stats_deriv_dir = options.stats_dir;
second_lvl_stats_dir = fullfile(stats_deriv_dir, '2ndlevel');
if ~exist(second_lvl_stats_dir)
    mkdir(second_lvl_stats_dir)
end

%fmrwhy_batch_specify2ndlevel(stats_dir, con_fns, params)


% -------
% For each TASK and RUN and ECHO
% -------
% Loop through sessions, tasks, runs, echoes.
for t = 1:numel(tasks)
    task = tasks{t};
    for r = 1:numel(runs)
        run = runs{r};
        for e = 1:numel(echoes)
            echo = echoes{e};

            for i = 1:3
                if ~(strcmp(task, 'emotion') && strcmp(run, '1')) && i>1
                    continue;
                end

                stats_dir = fullfile(second_lvl_stats_dir, ['task-' task '_run-' run '_echo-' echo '_con-' num2str(i)]);
                if ~exist(stats_dir)
                    mkdir(stats_dir)
                end

                con_fns = {};
                for s = 1:numel(subs)
                    sub = subs{s};
                    run_dir_stats = fullfile(stats_deriv_dir, ['sub-' sub], ['task-' task '_run-' run '_echo-' echo]);
                    fnc = fullfile(run_dir_stats, ['con_' sprintf('%04d', i) '_MNI152.nii']);
                    con_fns = [con_fns, {fnc}];
                end

                fmrwhy_batch_specify2ndlevel(stats_dir, con_fns, []);
                fmrwhy_batch_estimate1stlevel(stats_dir);
                consess = {};
                consess{1}.tcon.name = 'FingerTapping';
                consess{1}.tcon.weights = [1];
                consess{1}.tcon.sessrep = 'none';
                fmrwhy_batch_contrast1stlevel(stats_dir, consess)
                conspec = struct;
                conspec(1).titlestr = 'FingerTapping';
                conspec(1).contrasts = 1;
                conspec(1).threshdesc = 'none';
                conspec(1).thresh = 0.001;
                conspec(1).extent = 20;
                conspec(1).conjunction = 1;
                conspec(1).mask.none = 1;
                fmrwhy_batch_threshold2ndlevel(stats_dir, conspec)
            end
        end
    end
end

