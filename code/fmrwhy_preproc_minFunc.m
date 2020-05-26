function fmrwhy_preproc_minFunc(bids_dir, sub, ses, task, run, options)
%--------------------------------------------------------------------------

% Copyright statement....

%--------------------------------------------------------------------------
% DEFINITION
%--------------------------------------------------------------------------
% Function to run minimal functional preprocessing steps that are required for
% multi-echo combination

% STEPS:

% QUESTION: should functional localisers also be done based on combined echo data? Perhaps this is worth another research question?

% INPUT:

% OUTPUT:

%--------------------------------------------------------------------------

disp('---')
disp('*** Running fmrwhy_preproc_basicFunc ***')
disp('---')
disp('---')


% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

%% Update workflow params with subject functional derivative filenames
%options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

% Check if single / multi-echo (check from TE in settings; also possible to derive it from BIDS data)
is_multiecho = false;
if numel(options.TE) > 1
    is_multiecho = true;
end

% -------
% STEP 1: Estimate 3D volume realignment parameters from raw data
% -------
% First access template timeseries information
echo = options.template_echo;
options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);
stre_txt = ['sub-' sub '_task-' task '_run-' run '_echo-' echo];
% Check if realignment has already been done by seeing if the tsv file with head movement parameters exist
[d, f, e] = fileparts(options.motion_fn);
if ~exist(options.motion_fn, 'file')
    % If it does not exist estimate MPs
    disp(['Estimating 3D realignment parameters for: ' stre_txt]);
    realign_measures = fmrwhy_batch_realignEst(options.functional_fn, options.template_fn);
    temp_txt_fn = fullfile(d, [f '.txt']);
    col_names = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
    data = load(realign_measures.mp_fn);
    data_table = array2table(data,'VariableNames', col_names);
    writetable(data_table, temp_txt_fn, 'Delimiter','\t');
    [status, msg, msgID] = movefile(temp_txt_fn, options.motion_fn);
    disp('Complete!')
    disp('---')
else
    disp(['3D realignment parameters already estimated: ' options.motion_fn])
    disp('---')
end

%%
% -------
% STEP 2: Slice timing correction
% -------
if is_multiecho
    for e = 1:numel(options.TE)
        disp('---')
        disp(['Echo ' num2str(e)])
        disp('---')

        % Update workflow params with subject functional derivative filenames
        options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, num2str(e), options);
        % Filler text
        stre_txt = ['sub-' sub '_task-' task '_run-' run '_echo-' num2str(e)];
        if ~exist(options.afunctional_fn, 'file')
            disp(['Performing slice timing correction on: ' stre_txt])
            fmrwhy_batch_sliceTiming(options.functional_fn, options.afunctional_fn, options);
            disp('Complete!')
            disp('---')
        else
            disp(['Slice timing correction already completed for: ' stre_txt])
            disp('---')
        end
    end
else
    % TODO: implement automatic single-echo processing
end

%%
% -------
% STEP 3: 3D volume realignment
% -------
motion_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' options.template_echo '_desc-confounds_motion.tsv']);
motion_struct = tdfread(motion_fn)
motion_params = struct2array(motion_struct);
if is_multiecho
    for e = 1:numel(options.TE)
        disp('---')
        disp(['Echo ' num2str(e)])
        disp('---')
        % Filler text
        stre_txt = ['sub-' sub '_task-' task '_run-' run '_echo-' num2str(e)];
        % Update workflow params with subject functional derivative filenames
        options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, num2str(e), options);
        % Realign raw timeseries data
        if ~exist(options.rfunctional_fn, 'file')
            disp(['Performing 3D realignment on raw timeseries: ' stre_txt])
            fmrwhy_util_applyTransform(options.functional_fn, motion_params, options.template_fn, options.rfunctional_fn)
            disp('Complete!')
            disp('---')
        else
            disp(['3D realignment already completed for raw timeseries: ' stre_txt])
            disp('---')
        end
        % Realign slice time corrected timeseries data
        if ~exist(options.rafunctional_fn, 'file')
            disp(['Performing 3D realignment on slice time corrected timeseries: ' stre_txt])
            fmrwhy_util_applyTransform(options.afunctional_fn, motion_params, options.template_fn, options.rafunctional_fn)
            disp('Complete!')
            disp('---')
        else
            disp(['3D realignment already completed for slice time corrected timeseries: ' stre_txt])
            disp('---')
        end
    end
else
    % TODO: implement automatic single-echo processing
    % fmrwhy_batch_realignEstResl(options.functional_fn, options.template_fn, options.rfunctional_fn);
    % fmrwhy_batch_realignEstResl(options.afunctional_fn, options.template_fn, options.rafunctional_fn);
end
