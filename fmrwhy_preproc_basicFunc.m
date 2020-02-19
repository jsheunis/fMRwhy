function fmrwhy_preproc_basicFunc(bids_dir, sub, ses, task, run, echo)
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


% Load/create required defaults
% bids_dir = '/Volumes/Stephan_WD/NEUFEPME_data_BIDS';
% sub = '001';
% ses = '';
% task = 'motor'; % changed for fingertapping experiment. TODO: change back. and update functioning.
% run = '1';
% echo = '2';

disp('Loading BIDS directory and files')
template_task = 'motor'; % changed for fingertapping experiment. TODO: change back. and update functioning.
template_run = '1';
template_echo = '2';
defaults.TR = 2;
defaults.N_slices = 34;
fwhm = 7;
stre_txt = ['sub-' sub '_task-' task '_run-' run '_echo-' echo];
str_txt = ['sub-' sub '_task-' task '_run-' run];

% Directory and content setup
deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);
if ~exist(sub_dir_preproc, 'dir')
    mkdir(sub_dir_preproc)
    sub_dir_BIDS = fullfile(bids_dir, ['sub-' sub]);
    copyfile(sub_dir_BIDS, sub_dir_preproc)
end

% Grab functional timeseries filename,
functional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_bold.nii']);

% Grab template filename
template_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_space-individual_bold.nii']);

% Grab standard output filenames
motion_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_desc-confounds_motion.tsv']);
afunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-apreproc_bold.nii']);
rfunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rpreproc_bold.nii']);
rafunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_bold.nii']);
sfunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-spreproc_bold.nii']);
srfunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srpreproc_bold.nii']);
srafunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srapreproc_bold.nii']);
basic_func_out_fns = {motion_fn, afunctional_fn, rfunctional_fn, rafunctional_fn, sfunctional_fn, srfunctional_fn, srafunctional_fn};

disp('Complete!')
disp('---')


% -------
% STEP 1: Estimate 3D volume realignment parameters from raw data
% -------
% Check if this has already been done by seeing if the tsv file with head movement parameters exist
motion_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_desc-confounds_motion.tsv']);
[d, f, e] = fileparts(motion_fn);
if ~exist(motion_fn, 'file')
    % If it does not exist estimate MPs
    disp(['Estimating 3D realignment parameters for: ' stre_txt]);
    realign_measures = fmrwhy_batch_realignEst(functional_fn, template_fn);
    temp_txt_fn = fullfile(d, [f '.txt']);
    col_names = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
    data = load(realign_measures.mp_fn);
    data_table = array2table(data,'VariableNames', col_names);
    writetable(data_table, temp_txt_fn, 'Delimiter','\t');
    [status, msg, msgID] = movefile(temp_txt_fn, motion_fn);
    disp('Complete!')
    disp('---')
else
    disp(['3D realignment parameters already estimated: ' motion_fn])
    disp('---')
end

%%
% -------
% STEP 2: Slice timing correction
% TODO: implement automatic multi-echo processing
% -------
afunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-apreproc_bold.nii']);
if ~exist(afunctional_fn, 'file')
    disp(['Performing slice timing correction on: ' stre_txt])
    fmrwhy_batch_sliceTiming(functional_fn, afunctional_fn, defaults);
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
rfunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rpreproc_bold.nii']);
if ~exist(rfunctional_fn, 'file')
    disp(['Performing 3D realignment on raw timeseries: ' stre_txt])
    fmrwhy_batch_realignEstResl(functional_fn, template_fn, rfunctional_fn);
    disp('Complete!')
    disp('---')
else
    disp(['3D realignment already completed for raw timeseries: ' stre_txt])
    disp('---')
end
% Realign slice time corrected timeseries data
rafunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_bold.nii']);
if ~exist(rafunctional_fn, 'file')
    disp(['Performing 3D realignment on raw timeseries: ' stre_txt])
    fmrwhy_batch_realignEstResl(afunctional_fn, template_fn, rafunctional_fn);
    disp('Complete!')
    disp('---')
else
    disp(['3D realignment already completed for raw timeseries: ' stre_txt])
    disp('---')
end

%%
% -------
% STEP 4: spatial smoothing
% TODO: implement automatic multi-echo processing
% -------
% Smooth raw timeseries data
sfunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-spreproc_bold.nii']);
if ~exist(sfunctional_fn, 'file')
    disp(['Performing spatial smoothing on raw timeseries: ' stre_txt])
    fmrwhy_batch_smooth(functional_fn, sfunctional_fn, fwhm);
    disp('Complete!')
    disp('---')
else
    disp(['Spatial smoothing already completed for raw timeseries: ' stre_txt])
    disp('---')
end
% Smooth realigned timeseries data
srfunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srpreproc_bold.nii']);
if ~exist(srfunctional_fn, 'file')
    disp(['Performing spatial smoothing on realigned timeseries: ' stre_txt])
    fmrwhy_batch_smooth(rfunctional_fn, srfunctional_fn, fwhm);
    disp('Complete!')
    disp('---')
else
    disp(['Spatial smoothing already completed for realigned timeseries: ' stre_txt])
    disp('---')
end
% Smooth realigned and slice time corrected timeseries data
srafunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srapreproc_bold.nii']);
if ~exist(srafunctional_fn, 'file')
    disp(['Performing spatial smoothing on realigned and slice time corrected timeseries: ' stre_txt])
    fmrwhy_batch_smooth(rafunctional_fn, srafunctional_fn, fwhm);
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
fmrwhy_preproc_generateMultRegr(bids_dir, sub, ses, task, run, echo)