function options = fmrwhy_bids_getAnatDerivs(bids_dir, sub, options)

    if isempty(options)
        options = fmrwhy_bids_setupQcSubDirs(bids_dir, sub, options);
    end

    % T1w filename
    if isempty(options.anat_template_session)
        [filename, filepath] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'ext', '_T1w.nii');
    else
        [filename, filepath] = fmrwhy_bids_constructFilename('anat', 'sub', sub, 'ses', options.anat_template_session, 'ext', '_T1w.nii');
    end
    options.anatomical_fn = fullfile(options.preproc_dir, filepath, filename);
    % Outputs after coregistering T1w image
    options.coregest_anatomical_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'coregEst', 'ext', '_T1w.nii'));
    % Outputs after segmenting coregistered T1w image
    options.gm_probseg_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'GM', 'ext', '_probseg.nii'));
    options.wm_probseg_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'WM', 'ext', '_probseg.nii'));
    options.csf_probseg_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'CSF', 'ext', '_probseg.nii'));
    options.bone_probseg_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'bone', 'ext', '_probseg.nii'));
    options.soft_probseg_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'soft', 'ext', '_probseg.nii'));
    options.air_probseg_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'air', 'ext', '_probseg.nii'));
    options.indiv_to_mni_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'desc', 'IndivToMNI', 'ext', '_transform.nii')); % forward transform
    options.mni_to_indiv_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'desc', 'MNItoIndiv', 'ext', '_transform.nii')); % inverse transform
    options.probseg_fns = {options.gm_probseg_fn, options.wm_probseg_fn, options.csf_probseg_fn, options.bone_probseg_fn, options.soft_probseg_fn, options.air_probseg_fn};
    options.transform_fns = {options.indiv_to_mni_fn, options.mni_to_indiv_fn};
    % Outputs after reslicing segments and coregistered T1w image
    options.rcoregest_anatomical_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'coregEstResl_T1w', 'ext', '.nii'));
    options.rgm_probseg_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'rGM', 'ext', '_probseg.nii'));
    options.rwm_probseg_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'rWM', 'ext', '_probseg.nii'));
    options.rcsf_probseg_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'rCSF', 'ext', '_probseg.nii'));
    options.rbone_probseg_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'rbone', 'ext', '_probseg.nii'));
    options.rsoft_probseg_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'rsoft', 'ext', '_probseg.nii'));
    options.rair_probseg_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'rair', 'ext', '_probseg.nii'));
    options.rprobseg_fns = {options.rgm_probseg_fn, options.rwm_probseg_fn, options.rcsf_probseg_fn, options.rbone_probseg_fn, options.rsoft_probseg_fn, options.rair_probseg_fn};
    options.rall_fns = {options.rcoregest_anatomical_fn, options.rgm_probseg_fn, options.rwm_probseg_fn, options.rcsf_probseg_fn, options.rbone_probseg_fn, options.rsoft_probseg_fn, options.rair_probseg_fn};
    % Outputs after creating masks
    options.gm_mask_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'GM', 'ext', '_mask.nii'));
    options.wm_mask_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'WM', 'ext', '_mask.nii'));
    options.csf_mask_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'CSF', 'ext', '_mask.nii'));
    options.brain_mask_fn = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', 'brain', 'ext', '_mask.nii'));
    options.mask_fns = {options.gm_mask_fn, options.wm_mask_fn, options.csf_mask_fn, options.brain_mask_fn};
    % Outputs after anatomical localisation
    if options.map_rois == 1
        for i = 1:numel(options.tasks)
            % Ignore the 'rest' task (assume there is no task ROI for this; have to change in future if RSnetworks available to be normalised or something)
            if strcmp(options.tasks{i}, 'rest') ~= 1
                % Loop through all ROIs for the particular task
                for j = 1:numel(options.roi.(options.tasks{i}).orig_fn)
                    options.roi.(options.tasks{i}).roi_fn{j} = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', options.roi.(options.tasks{i}).desc{j}, 'ext', '_roi.nii'));
                    options.roi.(options.tasks{i}).rroi_fn{j} = fullfile(options.preproc_dir, filepath, fmrwhy_bids_constructFilename('anat', 'sub', sub, 'space', 'individual', 'desc', ['r' options.roi.(options.tasks{i}).desc{j}], 'ext', '_roi.nii'));
                end
            end
        end
    end
