function fmrwhy_preproc_anatLocaliser(bids_dir, sub, opts)
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


% Setup fmrwhy bids directories on workflow level
fmrwhy_defaults_setupDerivDirs(bids_dir);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
fmrwhy_defaults_setupSubDirs(bids_dir, sub);

% Update workflow params with subject anatomical derivative filenames
opts = fmrwhy_defaults_subAnat(bids_dir, sub, opts);


% -------
% STEP 1: Grab inputs necessary for normalisation (template and transform)
% -------
% Template functional volume
template_fn = fullfile(opts.sub_dir_preproc, 'func', ['sub-' sub '_task-' opts.template_task '_run-' opts.template_run '_space-individual_bold.nii']);
% Transformation from mni to individual functional template space
transformation_fn = opts.mni_to_indiv_fn;

%%
% -------
% STEP 2a: Write raw roi filenames to a cell array; construct normalised output filenames in BIDS format
% -------
toTransform_fns = {};
saveAs_fns = {};
count = 0;
% Loop through all tasks in BIDS structure
for i = 1:numel(opts.tasks)
    % Ignore the 'rest' task (assume there is no task ROI for this; have to change in future if RSnetworks available to be normalised or something)
    if strcmp(opts.tasks{i}, 'rest') ~= 1
        % Loop through all ROIs for the particular task
        for j = 1:numel(opts.roi.(opts.tasks{i}).orig_fn)
            count = count + 1;
            toTransform_fns{count} = opts.roi.(opts.tasks{i}).orig_fn{j};
            saveAs_fns{count} = fullfile(opts.anat_dir_preproc, ['sub-' sub '_space-individual_desc-' opts.roi.(opts.tasks{i}).desc{j} '_roi.nii']);
            opts.roi.(opts.tasks{i}).roi_fn{j} = saveAs_fns{count}; % save normalised filename for future use
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
for i = 1:numel(opts.tasks)
    % Ignore the 'rest' task (assume there is no task ROI for this; have to change in future if RSnetworks available to be normalised or something)
    if strcmp(opts.tasks{i}, 'rest') ~= 1
        % Loop through all ROIs for the particular task
        for j = 1:numel(opts.roi.(opts.tasks{i}).orig_fn)
            count = count + 1;
            reslice_fns{count} = opts.roi.(opts.tasks{i}).roi_fn{j};
            saveAs_fns{count} = fullfile(opts.anat_dir_preproc, ['sub-' sub '_space-individual_desc-r' opts.roi.(opts.tasks{i}).desc{j} '_roi.nii']);
            opts.roi.(opts.tasks{i}).rroi_fn{j} = saveAs_fns{count};% save resliced normalised filename for future use
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
% STEP 4: Show overlays to check quality
% TODO: should this be in a qc function/workflow?
% -------
count = 0;
[p1, frm1, rg1, dim1] = fmrwhy_util_readNifti(template_fn);
% Loop through all tasks in BIDS structure
for i = 1:numel(opts.tasks)
    % Ignore the 'rest' task (assume there is no task ROI for this; have to change in future if RSnetworks available to be normalised or something)
    if strcmp(opts.tasks{i}, 'rest') ~= 1
        % Loop through all ROIs for the particular task
        for j = 1:numel(opts.roi.(opts.tasks{i}).orig_fn)
            count = count+1;
            [p2, frm2, rg2, dim2] = fmrwhy_util_readNifti(opts.roi.(opts.tasks{i}).rroi_fn{j});
            overlay_img = fmrwhy_util_createBinaryImg(p2.nii.img, 0);
            title = opts.roi.(opts.tasks{i}).name{j}
            fmrwhy_util_createOverlayMontage(p1.nii.img, overlay_img, 9, 1, title, 'gray', 'on', 'max');
        end
    end
end