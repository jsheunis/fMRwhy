function defaults = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, defaults)

if isempty(defaults)
    defaults = fmrwhy_defaults_workflow(bids_dir)
end
% Subject level directories
defaults.sub_dir_preproc = fullfile(defaults.preproc_dir, ['sub-' sub]);
defaults.sub_dir_qc = fullfile(defaults.qc_dir, ['sub-' sub]);
defaults.func_dir_preproc = fullfile(defaults.sub_dir_preproc, 'func');
defaults.func_dir_qc = fullfile(defaults.sub_dir_qc, 'func');

% Template filename
defaults.template_fn = fullfile(defaults.sub_dir_preproc, 'func', ['sub-' sub '_task-' defaults.template_task '_run-' defaults.template_run '_space-individual_bold.nii']);
% Outputs from basicFunc processing
defaults.motion_fn = fullfile(defaults.func_dir_preproc, ['sub-' sub '_task-' defaults.template_task '_run-' defaults.template_run '_desc-confounds_motion.tsv']);
defaults.functional_fn = fullfile(defaults.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_bold.nii']);
defaults.afunctional_fn = fullfile(defaults.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-apreproc_bold.nii']);
defaults.rfunctional_fn = fullfile(defaults.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rpreproc_bold.nii']);
defaults.rafunctional_fn = fullfile(defaults.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_bold.nii']);
defaults.sfunctional_fn = fullfile(defaults.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-spreproc_bold.nii']);
defaults.srfunctional_fn = fullfile(defaults.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srpreproc_bold.nii']);
defaults.srafunctional_fn = fullfile(defaults.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srapreproc_bold.nii']);
defaults.framewise_displacement_fn = fullfile(defaults.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-confounds_fd.tsv']);
defaults.tissue_regr_fn = fullfile(defaults.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-confounds_tissue.tsv']);
defaults.physio_regr_fn = fullfile(defaults.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-confounds_physio.tsv']);
defaults.confounds_fn = fullfile(defaults.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-confounds_regressors.tsv']);
defaults.basic_func_out_fns = {defaults.motion_fn, defaults.afunctional_fn, defaults.rfunctional_fn, defaults.rafunctional_fn, defaults.sfunctional_fn, defaults.srfunctional_fn, defaults.srafunctional_fn, defaults.confounds_fn};