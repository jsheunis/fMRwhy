function fmrwhy_preproc_anatLocaliser(bids_dir, sub, options)
% FUNCTION:
%--------------------------------------------------------------------------

% Copyright statement....

%--------------------------------------------------------------------------
% DEFINITION
%--------------------------------------------------------------------------

%...

% INPUT:
% funcional0_fn     - filename of initial pre-real-time 3D functional volume template
% structural_fn     - filename of T1-weighted structural volume
% spm_dir           - SPM12 directory
%
% OUTPUT:
% output            - structure with filenames and data

%--------------------------------------------------------------------------
% STEPS

% Assumes basic structura-functional preproc has been done (coreg, etc), i.e.
% that the subject to MNI space mapping exists in the preproc_data stucture.
%--------------------------------------------------------------------------


disp('---')
disp('*** Running fmrwhy_preproc_anatLocaliser ***')
disp('---')
disp('---')

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);



% -------
% STEP 1: Grab inputs necessary for normalisation (template and transform)
% -------
% Template functional volume
template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
% Transformation from mni to individual functional template space
transformation_fn = options.mni_to_indiv_fn;

%%
% -------
% STEP 2a: Write raw roi filenames to a cell array; construct normalised output filenames in BIDS format
% -------
toTransform_fns = {};
saveAs_fns = {};
count = 0;
% Loop through all tasks in BIDS structure
for i = 1:numel(options.tasks)
    % Ignore the 'rest' task (assume there is no task ROI for this; have to change in future if RSnetworks available to be normalised or something)
    if strcmp(options.tasks{i}, 'rest') ~= 1
        % Loop through all ROIs for the particular task
        for j = 1:numel(options.roi.(options.tasks{i}).orig_fn)
            count = count + 1;
            toTransform_fns{count} = options.roi.(options.tasks{i}).orig_fn{j};
            saveAs_fns{count} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-' options.roi.(options.tasks{i}).desc{j} '_roi.nii']);
            options.roi.(options.tasks{i}).roi_fn{j} = saveAs_fns{count}; % save normalised filename for future use
        end
    end
end

%%
% -------
% STEP 2b:  Normalise all rois to individual functional template space
% -------
fmrwhy_batch_normaliseWrite(toTransform_fns, transformation_fn, template_fn, saveAs_fns)

%%
% -------
% STEP 3a:  Write raw roi filenames to a cell array; construct normalised output filenames in BIDS format
% -------
reslice_fns = {};
count = 0;
% Loop through all tasks in BIDS structure
for i = 1:numel(options.tasks)
    % Ignore the 'rest' task (assume there is no task ROI for this; have to change in future if RSnetworks available to be normalised or something)
    if strcmp(options.tasks{i}, 'rest') ~= 1
        % Loop through all ROIs for the particular task
        for j = 1:numel(options.roi.(options.tasks{i}).orig_fn)
            count = count + 1;
            reslice_fns{count} = options.roi.(options.tasks{i}).roi_fn{j};
            saveAs_fns{count} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-r' options.roi.(options.tasks{i}).desc{j} '_roi.nii']);
            options.roi.(options.tasks{i}).rroi_fn{j} = saveAs_fns{count};% save resliced normalised filename for future use
        end
    end
end

%%
% -------
% STEP 3b:  Reslice to functional-resolution image grid
% -------
fmrwhy_batch_coregResl(reslice_fns, template_fn, saveAs_fns)


%%
% -------
% STEP 4:  Create overlay montages
% -------
count = 0;
[p1, frm1, rg1, dim1] = fmrwhy_util_readOrientNifti(options.rcoregest_anatomical_fn);
% Loop through all tasks in BIDS structure
for i = 1:numel(options.tasks)
    % Ignore the 'rest' task (assume there is no task ROI for this; have to change in future if RSnetworks available to be normalised or something)
    if strcmp(options.tasks{i}, 'rest') ~= 1
        % Loop through all ROIs for the particular task
        for j = 1:numel(options.roi.(options.tasks{i}).orig_fn)
            count = count+1;
            [p2, frm2, rg2, dim2] = fmrwhy_util_readOrientNifti(options.roi.(options.tasks{i}).rroi_fn{j});
            overlay_img = fmrwhy_util_createBinaryImg(p2.nii.img, 0.1);
            title = options.roi.(options.tasks{i}).name{j}
            saveAs_fn = fullfile(options.anat_dir_qc, ['sub-' sub '_space-individual_desc-' options.roi.(options.tasks{i}).desc{j} '_roi_montage.png']);
            fmrwhy_util_createOverlayMontage(p1.nii.img, overlay_img, 9, 1, title, 'gray', 'off', 'max', [], [255, 0, 0], saveAs_fn);
        end
    end
end



disp('---')
disp('*** Finished running fmrwhy_preproc_anatLocaliser! ***')
disp('---')
disp('---')