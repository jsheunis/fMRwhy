% A custom workflow that does structFunc and basicFunc preprocessing and QC for a single subject in the NEUFEP study

% Code steps:
% 1. Define template/default variables, directories and filenames
% 2. Specify


%--------------------------------------------------------------------------


% -------
% STEP 0.1 -- Load defaults, filenames and parameters
% -------

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
% STEP 0.2 -- Create functional template
% -------
% Create, if it does not exist
template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
if ~exist(template_fn, 'file')
    disp(['Template funcional image does not exist yet. Creating now: ' template_fn]);
    functional0_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_echo-' options.template_echo '_bold.nii,1']);
    fmrwhy_util_saveNifti(template_fn, spm_read_vols(spm_vol(functional0_fn)), functional0_fn)
else
    disp(['Template funcional image exists: ' template_fn]);
end
options.template_fn = template_fn;



% -------
% PREPROCESSING PER TASK AND RUN
% -------

% Loop through sessions, tasks, runs, etc
ses = '';
tasks = {'rest', 'motor'};


for t = 1:numel(tasks)

    task = tasks{t};
    run = '1';
    echo = '2';

    % Update workflow params with subject functional derivative filenames
    options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

    % -------
    % STEP 1 -- Structural-functional preprocessing: fmrwhy_preproc_structFunc.m
    % -------
    % Loop through all standard output filenames and see if these files exist
    struct_func_out_fns = [{options.coregest_anatomical_fn} options.probseg_fns options.transform_fns options.rall_fns options.mask_fns];
    run_structFunc = 0;
    for i = 1:numel(struct_func_out_fns)
        if ~exist(struct_func_out_fns{i}, 'file')
            disp(['Structural-funcional preprocessing output file does not exist yet: ' struct_func_out_fns{i}]);
            run_structFunc = 1;
        end
    end
    % If some of the files do not exist, run the fmrwhy_preproc_structFunc processing pipeline
    if run_structFunc
        disp('Running complete structural-funcional preprocessing pipeline')
        fmrwhy_preproc_structFunc(bids_dir, sub, ses, options.template_task, options.template_run, options.template_echo, options);
        disp('Complete!')
        disp('---')
    else
        disp('Structural-funcional preprocessing already completed.')
        disp('---')
    end


    % -------
    % STEP 2 -- Basic functional preprocessing: fmrwhy_preproc_basicFunc.m
    % -------
    % Loop through all standard output filenames and see if these files exist
    basic_func_out_fns = {options.motion_fn, options.afunctional_fn, options.rfunctional_fn, options.rafunctional_fn, options.sfunctional_fn, options.srfunctional_fn, options.srafunctional_fn, options.confounds_fn};
    run_basicFunc = 0;
    for i = 1:numel(basic_func_out_fns)
        if ~exist(basic_func_out_fns{i}, 'file')
            disp(['Basic funcional preprocessing output file does not exist yet: ' basic_func_out_fns{i}]);
            run_basicFunc = 1;
        end
    end
    % If some of the files do not exist, run the fmrwhy_preproc_basicFunc processing pipeline
    if run_basicFunc
        fmrwhy_preproc_basicFunc(bids_dir, sub, ses, task, run, echo, options);
        disp('Complete!')
        disp('---')
    else
        disp('Basic funcional preprocessing already completed.')
        disp('---')
    end


    % -------
    % STEP 3 -- Anatomical localiser
    % TODO: add more checks to see if this was already done, and add logic to decide what to do
    % -------
    if options.map_rois == 1
        rrightAmygdala_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rrightAmygdala_roi.nii']);
        anatLocaliser_fns = {rrightAmygdala_fn};
        run_anatLocaliser = 0;
        for i = 1:numel(anatLocaliser_fns)
            if ~exist(anatLocaliser_fns{i}, 'file')
                disp(['Anatomical localiser output file does not exist yet: ' anatLocaliser_fns{i}]);
                run_anatLocaliser = 1;
            end
        end
        % If some of the files do not exist, run the fmrwhy_preproc_basicFunc processing pipeline
        if run_anatLocaliser
            fmrwhy_preproc_anatLocaliser(bids_dir, sub, options)
            disp('Complete!')
            disp('---')
        else
            disp('Anatomical localiser processing already completed.')
            disp('---')
        end

    end

    %%
    % -------
    % STEP 4 -- Quality control pipeline: fmrwhy_qc_run.m
    % -------

    fmrwhy_qc_run(bids_dir, sub, ses, task, run, echo, options);
    %qc_out_fns; % some file that is generated by the qc procedure (subject level or run level?)
    %run_qc = 0;
    %for i = 1:numel(qc_out_fns)
    %    if ~exist(qc_out_fns{i}, 'file')
    %        disp(['Basic funcional preprocessing output file does not exist yet: ' qc_out_fns{i}]);
    %        run_qc = 1;
    %    end
    %end
    %if run_qc
    %    fmrwhy_qc_run(bids_dir, sub, ses, task, run, echo, options);
    %    disp('Complete!')
    %    disp('---')
    %else
    %    disp('Basic funcional preprocessing already completed.')
    %    disp('---')
    %end


    %%
    % -------
    % STEP 5 -- QC report: fmrwhy_qc_generateSubRunReport.m
    % -------
    fmrwhy_qc_generateSubRunReport(bids_dir, sub, task, run, options)





    %
    %% Step 3: anatomical-localizer-preproc:     - rtme_preproc_anatLocaliser.m
    %fmrwhy_preproc_anatLocaliser(sub, options);

    % Step 3: functional-localizer-preproc:     - rtme_preproc_funcLocaliser.m
    %                                           - rtme_preproc_generateRegressors.m
    %                                           - rtme_preproc_generateRetroicor.m
    %                                           - rtme_preproc_generateFDregr.m
    %                                           - rtme_preproc_generateTissueSignals.m

    %
    %for t = 1:numel(defaults.tasks)
    %%    disp(['Performing 3D volume realignment for: ' sub '_task-' tasks(t) '_run-' defaults.template_run])
    %    rtme_preproc_funcLocaliser(sub, task, defaults.template_run, options.template_echo, defaults)
    %end


    % Step 4: calculate-prior-measures-preproc - rtme_preproc_estimateParams.m
    %rtme_preproc_estimateParams(sub, defaults);
end



