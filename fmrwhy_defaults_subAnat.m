function options = fmrwhy_defaults_subAnat(bids_dir, sub, options)

if isempty(options)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options)
end

% Subject level directories
options.anat_dir_preproc = fullfile(options.sub_dir_preproc, 'anat');
options.anat_dir_qc = fullfile(options.sub_dir_qc, 'anat');

% T1w filename
options.anatomical_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_T1w.nii']);
% Outputs after coregistering T1w image
options.coregest_anatomical_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-coregEst_T1w.nii']);
% Outputs after segmenting coregistered T1w image
options.gm_probseg_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-GM_probseg.nii']);
options.wm_probseg_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-WM_probseg.nii']);
options.csf_probseg_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-CSF_probseg.nii']);
options.bone_probseg_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-bone_probseg.nii']);
options.soft_probseg_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-soft_probseg.nii']);
options.air_probseg_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-air_probseg.nii']);
options.indiv_to_mni_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_desc-IndivToMNI_transform.nii']); % forward transform
options.mni_to_indiv_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_desc-MNItoIndiv_transform.nii']); % inverse transform
options.probseg_fns = {options.gm_probseg_fn, options.wm_probseg_fn, options.csf_probseg_fn, options.bone_probseg_fn, options.soft_probseg_fn, options.air_probseg_fn};
options.transform_fns = {options.indiv_to_mni_fn, options.mni_to_indiv_fn};
% Outputs after reslicing segments and coregistered T1w image
options.rcoregest_anatomical_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-coregEstResl_T1w.nii']);
options.rgm_probseg_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rGM_probseg.nii']);
options.rwm_probseg_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rWM_probseg.nii']);
options.rcsf_probseg_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rCSF_probseg.nii']);
options.rbone_probseg_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rbone_probseg.nii']);
options.rsoft_probseg_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rsoft_probseg.nii']);
options.rair_probseg_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rair_probseg.nii']);
options.rprobseg_fns = {options.rgm_probseg_fn, options.rwm_probseg_fn, options.rcsf_probseg_fn, options.rbone_probseg_fn, options.rsoft_probseg_fn, options.rair_probseg_fn};
options.rall_fns = {options.rcoregest_anatomical_fn, options.rgm_probseg_fn, options.rwm_probseg_fn, options.rcsf_probseg_fn, options.rbone_probseg_fn, options.rsoft_probseg_fn, options.rair_probseg_fn};
% Outputs after creating masks
options.gm_mask_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-GM_mask.nii']);
options.wm_mask_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-WM_mask.nii']);
options.csf_mask_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-CSF_mask.nii']);
options.brain_mask_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-brain_mask.nii']);
options.mask_fns = {options.gm_mask_fn, options.wm_mask_fn, options.csf_mask_fn, options.brain_mask_fn};