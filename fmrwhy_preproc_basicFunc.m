function fmrwhy_preproc_basicFunc(bids_dir, sub, ses, task, run, echo, options)
%--------------------------------------------------------------------------

% Copyright statement....

%--------------------------------------------------------------------------
% DEFINITION
%--------------------------------------------------------------------------
% Function to run basic functional preprocessing steps that are required for
% several subsequent analysis steps and quality control.

% THIS PIPELINE IS RUN ON A SIGNLE FUNCTIONAL TIMESERIES FILE
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

% Update workflow params with subject functional derivative filenames
options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

% Filler text
stre_txt = ['sub-' sub '_task-' task '_run-' run '_echo-' echo];
str_txt = ['sub-' sub '_task-' task '_run-' run];


% -------
% STEP 1: Estimate 3D volume realignment parameters from raw data
% TODO: implement automatic multi-echo processing
% -------
% Check if this has already been done by seeing if the tsv file with head movement parameters exist
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
% TODO: implement automatic multi-echo processing
% -------
if ~exist(options.afunctional_fn, 'file')
    disp(['Performing slice timing correction on: ' stre_txt])
    fmrwhy_batch_sliceTiming(options.functional_fn, options.afunctional_fn, options);
    disp('Complete!')
    disp('---')
else
    disp(['Slice timing correction already completed for: ' stre_txt])
    disp('---')
end

%%
% -------
% STEP 3: 3D volume realignment
% TODO: implement automatic multi-echo processing
% -------
% Realign raw timeseries data
if ~exist(options.rfunctional_fn, 'file')
    disp(['Performing 3D realignment on raw timeseries: ' stre_txt])
    fmrwhy_batch_realignEstResl(options.functional_fn, options.template_fn, options.rfunctional_fn);
    disp('Complete!')
    disp('---')
else
    disp(['3D realignment already completed for raw timeseries: ' stre_txt])
    disp('---')
end
% Realign slice time corrected timeseries data
if ~exist(options.rafunctional_fn, 'file')
    disp(['Performing 3D realignment on slice time corrected timeseries: ' stre_txt])
    fmrwhy_batch_realignEstResl(options.afunctional_fn, options.template_fn, options.rafunctional_fn);
    disp('Complete!')
    disp('---')
else
    disp(['3D realignment already completed for slice time corrected timeseries: ' stre_txt])
    disp('---')
end

%%
% -------
% STEP 4: spatial smoothing
% TODO: implement automatic multi-echo processing
% -------
% Smooth raw timeseries data
if ~exist(options.sfunctional_fn, 'file')
    disp(['Performing spatial smoothing on raw timeseries: ' stre_txt])
    fmrwhy_batch_smooth(options.functional_fn, options.sfunctional_fn, options.fwhm);
    disp('Complete!')
    disp('---')
else
    disp(['Spatial smoothing already completed for raw timeseries: ' stre_txt])
    disp('---')
end
% Smooth realigned timeseries data
if ~exist(options.srfunctional_fn, 'file')
    disp(['Performing spatial smoothing on realigned timeseries: ' stre_txt])
    fmrwhy_batch_smooth(options.rfunctional_fn, options.srfunctional_fn, options.fwhm);
    disp('Complete!')
    disp('---')
else
    disp(['Spatial smoothing already completed for realigned timeseries: ' stre_txt])
    disp('---')
end
% Smooth realigned and slice time corrected timeseries data
if ~exist(options.srafunctional_fn, 'file')
    disp(['Performing spatial smoothing on realigned and slice time corrected timeseries: ' stre_txt])
    fmrwhy_batch_smooth(options.rafunctional_fn, options.srafunctional_fn, options.fwhm);
    disp('Complete!')
    disp('---')
else
    disp(['Spatial smoothing already completed for realigned and slice time corrected timeseries: ' stre_txt])
    disp('---')
end


% -------
% STEP 5: Generate multiple regressors for GLM analysis and QC
% Includes: 3D realignment parameters, framewise displacement, FD censoring, tissue compartment signals, retroicor and HRV+RVT
% -------
if ~exist(options.confounds_fn, 'file')
    disp(['Generating multiple regressors for GLM analysis and QC'])
    fmrwhy_preproc_generateMultRegr(bids_dir, sub, ses, task, run, echo, options);
    disp('Complete!')
    disp('---')
else
    disp(['Multiple regressors already generated: ' options.confounds_fn])
    disp('---')
end
