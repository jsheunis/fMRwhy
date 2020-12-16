% A custom workflow that does ...

% Code steps:
% 1. Define template/default variables, directories and filenames
% 2. Create functional template, if it does not exist
% 3. Estimate 3D volume realignment parameters from raw template echo timeseries (given supplied template volume)
% 4. Run slice time correction for each echo timeseries
% 5. Realign each echo timeseries by applying rigid body transormation estimated from template echo realignment parameters
% 6. Calculate tSNR per echo, using the slice time corrected and realigned functional timeseries as inputs
% 7. Estimate T2star and S0 maps from minimally preprocessed multi-echo data

%--------------------------------------------------------------------------


% -------
% STEP 1 -- Load defaults, filenames and parameters
% -------
% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Set subject, sessions
sub = '001';
ses = '';

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);


% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
% Run template process on specific task and run predefined as template for multi-echo purposes
% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------

% -------
% STEP 5: Calculate/estimate T2star and S0 maps
% -------
% Template info
task = 'rest';
run = '1';

me_fns = {};
t2star_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_t2star.nii']);
s0_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_s0.nii']);

for e = 1:options.Ne
    % Update workflow params with subject functional derivative filenames
    options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, num2str(e), options);
    me_fns{e} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' num2str(e) '_desc-rapreproc_bold.nii']);
end
mask_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-brain_mask.nii']);
MEparams = fmrwhy_util_estimateMEparams(me_fns, options.TE, mask_fn, options.template_fn, t2star_fn, s0_fn);


overlaymontage = fmrwhy_util_createMontage(MEparams.T2star_3D_thresholded, 9, 1, 'xxx', 'hot', 'on', 'max', [0 200]); colorbar;
overlaymontage = fmrwhy_util_createMontage(MEparams.T2star_3D, 9, 1, 'xxx', 'hot', 'on', 'max', []); colorbar;


[p, frm, rg, dim] = fmrwhy_util_readOrientNifti(t2star_fn);
class(p.nii.img)
img = double(p.nii.img);
class(img)
overlaymontage = fmrwhy_util_createMontage(img, 9, 1, 'xxx', 'hot', 'on', 'max', []); colorbar;