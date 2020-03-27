% A custom workflow that does structFunc and basicFunc preprocessing and QC for a single subject in the NEUFEP study

% Code steps:
% 1. Define template/default variables, directories and filenames
% 2. Specify


%--------------------------------------------------------------------------


% -------
% STEP 1 -- Load defaults, filenames and parameters
% -------
disp('---')
disp('STEP 1 -- Load defaults, filenames and parameters')
disp('---')

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

% -------
% STEP 2 -- Create functional template, if it does not exist
% -------
disp('---')
disp('STEP 2 -- Create functional template, if it does not exist')
disp('---')
% Create, if it does not exist
template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
if ~exist(template_fn, 'file')
    disp(['Template funcional image does not exist yet. Creating now: ' template_fn]);
    functional0_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_echo-' options.template_echo '_bold.nii,1']);
    fmrwhy_util_saveNifti(template_fn, spm_read_vols(spm_vol(functional0_fn)), functional0_fn, 'Template functional volume', 0)
else
    disp(['Template funcional image exists: ' template_fn]);
end
options.template_fn = template_fn;


% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
% Run process on specific task and run predefined as template for multi-echo purposes
% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
ses = '';
task = 'rest';
run = '1';
% -------
% STEP 3: Estimate 3D volume realignment parameters from raw template echo timeseries
% -------
disp('---')
disp('STEP 3: Estimate 3D volume realignment parameters')
disp('---')
% Check if this has already been done by seeing if the tsv file with head movement parameters exist
echo = options.template_echo;
options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);
[d, f, e] = fileparts(options.motion_fn);
if ~exist(options.motion_fn, 'file')
    % If it does not exist estimate MPs
    disp('---')
    disp(['Estimating 3D realignment parameters for template echo timeseries']);
    realign_measures = fmrwhy_batch_realignEst(options.functional_fn, options.template_fn);
    temp_txt_fn = fullfile(d, [f '.txt']);
    col_names = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
    data = load(realign_measures.mp_fn);
    motion_params = data;
    data_table = array2table(data,'VariableNames', col_names);
    writetable(data_table, temp_txt_fn, 'Delimiter','\t');
    [status, msg, msgID] = movefile(temp_txt_fn, options.motion_fn);
    disp('Complete!')
else
    disp('---')
    disp(['3D realignment parameters already estimated, loading now: ' options.motion_fn])
    motion_struct = tdfread(options.motion_fn)
    motion_params = struct2array(motion_struct);
end
% -------
% STEP 4: slice time correction for each echo timeseries
% -------
disp('---')
disp('STEP 4: Slice time correction for each echo timeseries')
disp('---')
for e = 1:options.Ne
    echo = num2str(e);
    disp('---')
    disp(['Echo ' echo])
    disp('---')
    % Update workflow params with subject functional derivative filenames
    options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

    if ~exist(options.afunctional_fn, 'file')
        disp('---')
        disp(['Slice time corrected file does not exist yet: ' options.afunctional_fn]);
        fmrwhy_batch_sliceTiming(options.functional_fn, options.afunctional_fn, options)
        disp('Complete!')
        disp('---')
    else
        disp(['Slice time correction already completed: ' options.afunctional_fn])
        disp('---')
    end
end
% -------
% STEP 5: Realignment for each echo timeseries
% -------
disp('---')
disp('STEP 5: Realignment for each echo timeseries')
disp('---')
% For each echo apply transormation derived from motion parameters
for e = 1:options.Ne
    echo = num2str(e);
    disp('---')
    disp(['Echo ' echo])
    disp('---')
    % Update workflow params with subject functional derivative filenames
    options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

    if ~exist(options.rafunctional_fn, 'file')
        disp('---')
        disp(['Realigned + slice time corrected file does not exist yet: ' options.rafunctional_fn]);
        fmrwhy_util_applyTransform(options.afunctional_fn, motion_params, options.template_fn, options.rafunctional_fn)
        disp('Complete!')
        disp('---')
    else
        disp(['Realignment + slice time correction already completed: ' options.rafunctional_fn])
        disp('---')
    end
end


% -------
% STEP 6: Calculate tSNR per echo, using the slice time corrected and realigned functional timeseries!!!
% -------
disp('---')
disp('STEP 6: Calculate tSNR for each echo timeseries')
disp('---')
for e = 1:options.Ne
    echo = num2str(e);
    disp('---')
    disp(['Echo ' echo])
    disp('---')
    % Update workflow params with subject functional derivative filenames
    options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);


    mean_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_mean.nii']);
    std_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_std.nii']);
    tsnr_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_tsnr.nii']);

    if ~exist(options.rafunctional_fn, 'file')
        disp('---')
        disp(['Realigned + slice time corrected file does not exist yet: ' options.rafunctional_fn]);
        fmrwhy_util_applyTransform(options.afunctional_fn, motion_params, options.template_fn, options.rafunctional_fn)
        disp('Complete!')
        disp('---')
    else
        disp(['Realignment + slice time correction already completed: ' options.rafunctional_fn])
        disp('---')
    end
end


for e = 1:N_e
    disp(['tSNR for echo ' num2str(e)])
    F{e} = spm_read_vols(spm_vol(rf_me_fn{e}));
    F_ave2D{e} = mean(reshape(F{e},Ni*Nj*Nk, N_vol), 2);
    F_ave{e} = reshape(F_ave2D{e}, Ni, Nj, Nk);
    F_tSNR2D{e} = F_ave2D{e}./std(reshape(F{e},Ni*Nj*Nk, N_vol), 0, 2);
    F_tSNR{e} = reshape(F_tSNR2D{e}, Ni, Nj, Nk);
end



