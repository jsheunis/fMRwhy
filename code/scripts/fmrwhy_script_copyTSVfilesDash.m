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
bids_dir = '/Volumes/TSM/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Loop through subjects, sessions, tasks, runs, etc
%subs = {'016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
subs = {'001', '002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
%subs = {'001'};
ses = '';
%tasks = {'motor', 'emotion'};
%runs = {'1', '2'};
%echoes = {'2', 'combinedMEtsnr', 'combinedMEt2star', 'combinedMEte'};
tasks = {'motor', 'emotion'};
runs = {'1', '2'};
echoes = {'2', 'combinedMEtsnr', 'combinedMEt2star', 'combinedMEte'};
report_dir = '/Users/jheunis/Documents/Websites/rt-me-fmri-dash/assets/';

for s = 1:numel(subs)

    sub = subs{s};
    disp(sub);
    % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

    func_dir_me = fullfile(options.me_dir, ['sub-' sub], 'func');

    cd(func_dir_me);
    all_tsvs = dir('*.tsv');
    for i = 1:numel(all_tsvs)
        src = fullfile(func_dir_me, all_tsvs(i).name);
        copyfile(src, report_dir);
    end

%    transformation_fn = options.indiv_to_mni_fn;
%    template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
%    fmrwhy_batch_normaliseWrite(toTransform_fns, transformation_fn, template_fn, saveAs_transform_fns)
end