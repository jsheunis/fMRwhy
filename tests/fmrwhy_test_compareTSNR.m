bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';
sub = '001';
func_dir_preproc = fullfile(bids_dir, 'derivatives', 'fmrwhy-preproc', ['sub-' sub], 'func');
anat_dir_preproc = fullfile(bids_dir, 'derivatives', 'fmrwhy-preproc', ['sub-' sub], 'anat');

% Template details
task = 'rest';
run = '1';
template_fn = fullfile(func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_space-individual_bold.nii']);

% Mask details
masks = fmrwhy_util_loadMasks(bids_dir, sub);
mask_fn = masks.brain_mask_fn;

% Functional details
task = 'motor';
%task = 'emotion';
run = '1';
echo = '2';

% tSNR filenames
tsnr_fns = {};
tsnr_fns{1} = fullfile(func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_tsnr.nii']);
tsnr_fns{2} = fullfile(func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEt2star_tsnr.nii']);
tsnr_fns{3} = fullfile(func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEtsnr_tsnr.nii']);
tsnr_fns{4} = fullfile(func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEte_tsnr.nii']);

% ROIs
roi_fns = {};
roi_fns{1} = fullfile(anat_dir_preproc, ['sub-' sub '_space-individual_desc-rleftMotor_roi.nii']);
%roi_fns{1} = fullfile(anat_dir_preproc, ['sub-' sub '_space-individual_desc-rrightMotor_roi.nii']);
roi_text  = 'left motor cortex';
%roi_fns{1} = fullfile(anat_dir_preproc, ['sub-' sub '_space-individual_desc-rleftAmygdala_roi.nii']);
%roi_fns{2} = fullfile(anat_dir_preproc, ['sub-' sub '_space-individual_desc-rrightAmygdala_roi.nii']);
%roi_text  = 'bilateral amygdala';

fmrwhy_util_compareTSNR(tsnr_fns, mask_fn, roi_fns, roi_text)