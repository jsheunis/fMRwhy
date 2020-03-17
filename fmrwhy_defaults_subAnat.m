function defaults = fmrwhy_defaults_subAnat(bids_dir, sub, defaults)

if isempty(defaults)
    defaults = fmrwhy_defaults_workflow(bids_dir)
end
% Subject level directories
defaults.sub_dir_preproc = fullfile(defaults.preproc_dir, ['sub-' sub]);
defaults.sub_dir_qc = fullfile(defaults.qc_dir, ['sub-' sub]);
defaults.anat_dir_preproc = fullfile(defaults.sub_dir_preproc, 'anat');
defaults.anat_dir_qc = fullfile(defaults.sub_dir_qc, 'anat');
% T1w filename
defaults.anatomical_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_T1w.nii']);
% Outputs after coregistering T1w image
defaults.coregest_anatomical_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-coregEst_T1w.nii']);
% Outputs after segmenting coregistered T1w image
defaults.gm_probseg_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-GM_probseg.nii']);
defaults.wm_probseg_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-WM_probseg.nii']);
defaults.csf_probseg_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-CSF_probseg.nii']);
defaults.bone_probseg_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-bone_probseg.nii']);
defaults.soft_probseg_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-soft_probseg.nii']);
defaults.air_probseg_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-air_probseg.nii']);
defaults.indiv_to_mni_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_desc-IndivToMNI_transform.nii']); % forward transform
defaults.mni_to_indiv_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_desc-MNItoIndiv_transform.nii']); % inverse transform
defaults.probseg_fns = {defaults.gm_probseg_fn, defaults.wm_probseg_fn, defaults.csf_probseg_fn, defaults.bone_probseg_fn, defaults.soft_probseg_fn, defaults.air_probseg_fn};
defaults.transform_fns = {defaults.indiv_to_mni_fn, defaults.mni_to_indiv_fn};
% Outputs after reslicing segments and coregistered T1w image
defaults.rcoregest_anatomical_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-coregEstResl_T1w.nii']);
defaults.rgm_probseg_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rGM_probseg.nii']);
defaults.rwm_probseg_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rWM_probseg.nii']);
defaults.rcsf_probseg_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rCSF_probseg.nii']);
defaults.rbone_probseg_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rbone_probseg.nii']);
defaults.rsoft_probseg_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rsoft_probseg.nii']);
defaults.rair_probseg_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rair_probseg.nii']);
defaults.rprobseg_fns = {defaults.rgm_probseg_fn, defaults.rwm_probseg_fn, defaults.rcsf_probseg_fn, defaults.rbone_probseg_fn, defaults.rsoft_probseg_fn, defaults.rair_probseg_fn};
defaults.rall_fns = {defaults.rcoregest_anatomical_fn, defaults.rgm_probseg_fn, defaults.rwm_probseg_fn, defaults.rcsf_probseg_fn, defaults.rbone_probseg_fn, defaults.rsoft_probseg_fn, defaults.rair_probseg_fn};
% Outputs after creating masks
defaults.gm_mask_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-GM_mask.nii']);
defaults.wm_mask_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-WM_mask.nii']);
defaults.csf_mask_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-CSF_mask.nii']);
defaults.brain_mask_fn = fullfile(defaults.anat_dir_preproc, ['sub-' sub '_space-individual_desc-brain_mask.nii']);
defaults.mask_fns = {defaults.gm_mask_fn, defaults.wm_mask_fn, defaults.csf_mask_fn, defaults.brain_mask_fn};