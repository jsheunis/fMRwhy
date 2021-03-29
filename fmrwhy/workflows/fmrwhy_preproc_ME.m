function fmrwhy_preproc_ME(bids_dir, sub, ses, task, run, options)

% A workflow that does the minimal preprocessing necessary for multi-echo combination

% Code steps:
% 3. Estimate 3D volume realignment parameters from raw template echo timeseries (given supplied template volume)
% 4. Run slice time correction for each echo timeseries
% 5. Realign each echo timeseries by applying rigid body transormation estimated from template echo realignment parameters

%--------------------------------------------------------------------------



disp('---')
disp('*** Running fmrwhy_preproc_multiecho ***')
disp('---')
disp('---')


% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% -------
% STEP 1: Estimate 3D volume realignment parameters from raw timeseries (using given task, run, template echo)
% -------
disp('---')
disp('STEP 1: Estimate 3D volume realignment parameters')
disp('---')
% Check if this has already been done by seeing if the tsv file with head movement parameters exist
% First access template timeseries informatopn
options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, options.template_echo, options);
%options.motion_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' options.template_echo '_desc-confounds_motion.tsv'])
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
    disp('---')
    disp(['Echo ' num2str(e)])
    disp('---')
    % Update workflow params with subject functional derivative filenames
    options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, num2str(e), options);

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
% -------
disp('---')
disp('STEP 5: Realignment for each echo timeseries')
disp('---')
% For each echo apply transormation derived from motion parameters
for e = 1:options.Ne
    disp('---')
    disp(['Echo ' num2str(e)])
    disp('---')
    % Update workflow params with subject functional derivative filenames
    options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, num2str(e), options);

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


