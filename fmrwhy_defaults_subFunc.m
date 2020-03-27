function options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options)

if isempty(options)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);
end

% Template filename
options.template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
% Outputs from basicFunc processing
options.motion_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_echo-' options.template_echo '_desc-confounds_motion.tsv']);
options.functional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_bold.nii']);
options.afunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-apreproc_bold.nii']);
options.rfunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rpreproc_bold.nii']);
options.rafunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_bold.nii']);
options.sfunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-spreproc_bold.nii']);
options.srfunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srpreproc_bold.nii']);
options.srafunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srapreproc_bold.nii']);
options.framewise_displacement_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-confounds_fd.tsv']);
options.tissue_regr_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-confounds_tissue.tsv']);
options.physio_regr_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-confounds_physio.tsv']);
options.confounds_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-confounds_regressors.tsv']);
options.basic_func_out_fns = {options.motion_fn, options.afunctional_fn, options.rfunctional_fn, options.rafunctional_fn, options.sfunctional_fn, options.srfunctional_fn, options.srafunctional_fn, options.confounds_fn};