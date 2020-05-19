function fmrwhy_workflow_1stlevel(bids_dir, sub, ses, task, run, echo, options)

% A custom workflow that runs 1st level analysis for a single run based on specified parameters

% STEPS:

%--------------------------------------------------------------------------


% Load/create required defaults
% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% Update workflow params with subject functional derivative filenames
options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);


% -------
% STEP 1: set up directories for outputs
% -------
options.sub_dir_stats = fullfile(options.stats_dir, ['sub-' sub]);
func_dir_stats = fullfile(options.sub_dir_stats, 'func');
if ~exist(func_dir_stats, 'dir')
    mkdir(func_dir_stats)
end

% -------
% STEP 2: create design matrix regressors
% -------
% Load multiple confound regressors
confounds_struct = tdfread(options.confounds_fn);
confounds_mat = struct2array(confounds_struct);
regressors_mat = [];
regressors_names = {};

fields = fieldnames(options.firstlevel.glm_regressors);

for i = 1:numel(fields)
    key = fields{i};
    val = options.firstlevel.glm_regressors.(key);
    % If the value is true or larger than zero (for retroicor order), parse key and include relevant data in regressor matrix
    if val
        if strfind(key,'trans_rot')
            if strcmp(key,'trans_rot')
                trans_rot_keys = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
            elseif strcmp(key,'trans_rot_derivative1')
                trans_rot_keys = {'trans_x_derivative1', 'trans_y_derivative1', 'trans_z_derivative1', 'rot_x_derivative1', 'rot_y_derivative1', 'rot_z_derivative1'};
            elseif strcmp(key,'trans_rot_power2')
                trans_rot_keys = {'trans_x_power2', 'trans_y_power2', 'trans_z_power2', 'rot_x_power2', 'rot_y_power2', 'rot_z_power2'};
            elseif strcmp(key,'trans_rot_derivative1_power2')
                trans_rot_keys = {'trans_x_derivative1_power2', 'trans_y_derivative1_power2', 'trans_z_derivative1_power2', 'rot_x_derivative1_power2', 'rot_y_derivative1_power2', 'rot_z_derivative1_power2'};
            else
            end
            for j = 1:numel(trans_rot_keys)
                regressors_mat = [regressors_mat confounds_struct.(trans_rot_keys{j})];
                regressors_names = [regressors_names {trans_rot_keys{j}}];
            end
        elseif strcmp(key,'retroicor_c') || strcmp(key,'retroicor_r') || strcmp(key,'retroicor_cxr')
            for k=1:val
                new_key = [key num2str(k)];
                regressors_mat = [regressors_mat confounds_struct.(new_key)];
                regressors_names = [regressors_names {new_key}];
            end
        else
            regressors_mat = [regressors_mat confounds_struct.(key)];
            regressors_names = [regressors_names {key}];
        end
    end
end

regressors_fn = fullfile(func_dir_stats, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-GLM_regressors.txt']);
disp(regressors_names)
dlmwrite(regressors_fn, regressors_mat, 'delimiter', '\t', 'precision', '%1.7e')


% -------
% STEP 3: Set up statistical design parameters, based on task data
% -------
%% CREATE MODEL

events_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_events.tsv']);
events_struct = tdfread(events_fn)
%assignin('base','events_struct',events_struct)


options.firstlevel.(task).sess_params.timing_units = 'secs';
options.firstlevel.(task).sess_params.timing_RT = 2;
cond_names = {'Shapes', 'Faces'};
[cond, trials] = fmrwhy_util_1stlevelBIDStoConditions(events_fn, cond_names);


options.firstlevel.(task).sess_params.cond = cond;
sess_params = options.firstlevel.(task).sess_params;
%functional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srapreproc_bold.nii']);
functional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-scombinedMEt2star_bold.nii'])
fmrwhy_batch_specify1stlevel(func_dir_stats, functional_fn, regressors_fn, sess_params)
load([func_dir_stats filesep 'SPM.mat']);

%% ESTIMATE MODEL
fmrwhy_batch_estimate1stlevel(func_dir_stats)

%% SETUP TASK CONTRAST
[Ntt, Nregr] = size(SPM.xX.X);
contrast_params = struct;
contrast_params.weights = zeros(1, Nregr);
contrast_params.weights(2) = 1;
contrast_params.weights(1) = -1;
contrast_params.name = 'Faces';
fmrwhy_batch_contrast1stlevel(func_dir_stats, contrast_params)


% RUN RESULTS
fmrwhy_batch_threshold1stlevel(func_dir_stats)
[SPM, xSPM] = spm_getSPM(fullfile(func_dir_stats, 'SPM.mat'));
assignin('base', 'SPM', SPM)
