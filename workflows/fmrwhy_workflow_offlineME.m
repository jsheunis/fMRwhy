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

% Set subject, sessions
sub = '001';
ses = '';

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

% -------
% STEP 3: For all tasks and runs, complete minimal multi-echo preprocessing
% -------
disp('---')
disp('STEP 3: For all tasks and runs, complete minimal multi-echo preprocessing')
disp('---')
tasks = {'rest', 'motor', 'emotion'};
runs = {'1', '2'};

for t = 1:numel(tasks)
    task = tasks{t};
    for r = 1:numel(runs)
        run = runs{r};
        disp('------------')
        disp('------------')
        disp(['Task: ' task ';  Run: ' run])
        disp('------------')
        disp('------------')
        fmrwhy_preproc_ME(bids_dir, sub, ses, task, run, options)
    end
end

% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------
% Run template process on specific task and run predefined as template for multi-echo purposes
% -----------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------

% Template info
task = 'rest';
run = '1';

% -------
% STEP 4: Calculate tSNR per echo, using the slice time corrected and realigned functional timeseries!!!
% -------
disp('---')
disp('STEP 4: Calculate tSNR for each echo timeseries of template task+run')
disp('---')
for e = 1:options.Ne
    disp('---')
    disp(['Echo ' num2str(e)])
    disp('---')
    % Update workflow params with subject functional derivative filenames
    options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, num2str(e), options);
    tsnr_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' num2str(e) '_desc-rapreproc_tsnr.nii']);

    if ~exist(tsnr_fn, 'file')
        disp('---')
        disp(['Temporal SNR file does not exist yet. Calculating now...']);
        tsnr_output = fmrwhy_util_calculateTSNR(options.rafunctional_fn, 0, tsnr_fn, options.template_fn)
        disp('Complete!')
        disp('---')
    else
        disp(['Temporal SNR already calculated. Loading now...'])
        tsnr_output = struct;
        tsnr_output.data_3D_tsnr = spm_read_vols(spm_vol(tsnr_fn));
        disp('---')
    end
end


% -------
% STEP 5: Calculate/estimate T2star and S0 maps
% -------
disp('---')
disp('STEP 5: Calculate T2star and S0 maps')
disp('---')

t2star_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_t2star.nii']);
s0_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_s0.nii']);
if ~exist(t2star_fn, 'file')
    disp('---')
    disp(['T2star map file does not exist yet. Calculating now...']);
    me_fns = {};
    for e = 1:options.Ne
        % Update workflow params with subject functional derivative filenames
        options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, num2str(e), options);
        me_fns{e} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' num2str(e) '_desc-rapreproc_bold.nii']);
    end
    mask_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-brain_mask.nii']);
    MEparams = fmrwhy_util_estimateMEparams(me_fns, options.TE, mask_fn, options.template_fn, t2star_fn, s0_fn);
    disp('Complete!')
    disp('---')
else
    disp(['T2star map already calculated. Loading now...'])
    MEparams = struct;
    MEparams.T2star_3D_thresholded = spm_read_vols(spm_vol(t2star_fn));
    MEparams.S0_3D_thresholded = spm_read_vols(spm_vol(s0_fn));
    disp('---')
end

%% Visualise t2star and s0 maps
%[p1, frm1, rg1, dim1] = fmrwhy_util_readNifti(t2star_fn);
%[p2, frm2, rg2, dim2] = fmrwhy_util_readNifti(s0_fn);
%t2star_montage = fmrwhy_util_createMontage(p1.nii.img, 9, 1, 'T2star', 'hot', 'on', 'max');
%colorbar; % caxis([0 200]);
%s0_montage = fmrwhy_util_createMontage(p2.nii.img, 9, 1, 'S0', 'parula', 'on', 'max');
%colorbar;