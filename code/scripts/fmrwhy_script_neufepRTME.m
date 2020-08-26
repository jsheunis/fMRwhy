% ---------------------- %
% Main NEUFEP-ME real-time processing script
% ---------------------- %

% TODO: NOTE - EVERYTHING IS READ IN AND SAVED WITH SPM_VOL (ETC) AND NOT USING NII_TOOL


% -------
% STEP 1: Define run details
% -------
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

tasks = {'motor', 'emotion'};
tasks = {'emotion'};
runs = {'1', '2'};
runs = {'2'};
echoes = {'2', 'combinedMEtsnr', 'combinedMEt2star', 'combinedMEte'};
%subs = {'001', '002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
subs = {'001'};
ses = '';
echo = echoes{1};

% -------
% STEP 2: fMRwhy setup
% -------


for s = 1:numel(subs)

    clearvars -except bids_dir tasks runs echoes subs ses echo s

    sub = subs{s};

    options = struct;
    % Setup fmrwhy BIDS-derivatuve directories on workflow level
    options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);
    % Grab parameters from workflow settings file
    options = fmrwhy_settings_preprocQC(bids_dir, options);
    % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);
    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_defaults_subAnat(bids_dir, sub, options);
    % Multi-echo derivatives
    options.me_dir = fullfile(options.deriv_dir, 'fmrwhy-multiecho');
    options.sub_dir_me = fullfile(options.me_dir, ['sub-' sub]);
    % Create derivatives directory for rt output
    options.rt_dir = fullfile(options.deriv_dir, 'fmrwhy-rt');
    if ~exist(options.rt_dir, 'dir')
        mkdir(options.rt_dir);
    end
    % Create sub directory for rt output
    options.sub_dir_rt = fullfile(options.rt_dir, ['sub-' sub]);
    if ~exist(options.sub_dir_rt, 'dir')
        mkdir(options.sub_dir_rt);
    end

    for t = 1:numel(tasks)
        task = tasks{t};

        for r = 1:numel(runs)
            run = runs{r};

            disp('-----------------------------------------------------')
            disp(['Running RT analysis for: sub-' sub '_task-' task '_run-' run])
            disp('-----------------------------------------------------')

            % Update workflow params with subject functional derivative filenames
            options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);
            disp('rtme-init')
            fmrwhy_script_rtme_init;
            disp('rtme-script')
            fmrwhy_script_rtme;
            disp('rtme-postprocess')
            fmrwhy_script_rtme_postprocess;
        end
    end
end
