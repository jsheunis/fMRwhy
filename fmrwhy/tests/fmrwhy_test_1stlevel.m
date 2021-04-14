% A custom workflow that does ...

%% --------------------------------------------------------------------------
% clear;
%
%% Load/create required parameters
% bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';
%
%% Loop through subjects, sessions, tasks, runs, etc
% sub = '001';
% ses = '';
% task = 'motor';
%% tasks = {'motor', 'emotion'};
%% task = 'emotion';
% run = '1';
%% runs = {'1', '2'};
% echo = 'combinedMEt2star';
%% echoes = {'combinedMEt2star', 'combinedMEte'};
%
%% task = 'emotion'
%% run = '2'
%%% echo = 'combinedMEt2star'
%% echo = 'combinedMEte'
%
% options = struct;
% fmrwhy_workflow_1stlevelRun(bids_dir, sub, ses, task, run, echo, options)
%%
%% for t = 1:2
%%    task = tasks{t}
%%    for r = 1:2
%%        run = runs{t}
%%        for e = 1:2
%%            echo = echoes{e}
%%            options = struct;
%%            fmrwhy_workflow_1stlevelRun(bids_dir, sub, ses, task, run, echo, options)
%%        end
%%    end
%% end

bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';
sub = '001';
ses = '';
tasks = {'motor', 'emotion'};
runs = {'1', '2'};
echoes = {'2', 'combinedMEtsnr', 'combinedMEt2star', 'combinedMEte'};
options = struct;
% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);
% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);
% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);
% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

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

            for i = 1:3
                fn1 = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', i) '.nii']);
                fn2 = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', i) '_binary_clusters.nii']);
                fn3 = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', i) '_nary_clusters.nii']);

                fnsave1 = strrep(fn1, '.nii', '_MNI152.nii');
                fnsave2 = strrep(fn2, '.nii', '_MNI152.nii');
                fnsave3 = strrep(fn3, '.nii', '_MNI152.nii');

                if exist(fn1, 'file')
                    toTransform_fns = [toTransform_fns, {fn1, fn2, fn3}];
                    saveAs_transform_fns = [saveAs_transform_fns, {fnsave1, fnsave2, fnsave3}];
                else
                    continue
                end
            end

        end
    end
end

transformation_fn = options.indiv_to_mni_fn;
template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
fmrwhy_batch_normaliseWrite(toTransform_fns, transformation_fn, template_fn, saveAs_transform_fns);
