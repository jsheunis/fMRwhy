bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';
sub = '001';
func_dir_preproc = fullfile(bids_dir, 'derivatives', 'fmrwhy-preproc', ['sub-' sub], 'func');
anat_dir_preproc = fullfile(bids_dir, 'derivatives', 'fmrwhy-preproc', ['sub-' sub], 'anat');

% Template details
task = 'rest';
run = '1';
template_fn = fullfile(func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_space-individual_bold.nii']);

% Mask details
%masks = fmrwhy_util_loadMasks(bids_dir, sub);
%mask_fn = masks.brain_mask_fn;

% Functional details
task = 'motor';
%task = 'emotion';
run = '1';
%echo = '2';

% tSNR filenames
tsnr_fns = {};
tsnr_fns{1} = fullfile(func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_tsnr.nii']);
tsnr_fns{2} = fullfile(func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEt2star_tsnr.nii']);
tsnr_fns{3} = fullfile(func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEtsnr_tsnr.nii']);
tsnr_fns{4} = fullfile(func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEte_tsnr.nii']);

% ROIs
roi_fns = {};
%roi_fns{1} = fullfile(anat_dir_preproc, ['sub-' sub '_space-individual_desc-rleftMotor_roi.nii']);
roi_fns{1} = fullfile(anat_dir_preproc, ['sub-' sub '_space-individual_desc-rrightMotor_roi.nii']);
roi_text  = 'left motor cortex';
%roi_fns{1} = fullfile(anat_dir_preproc, ['sub-' sub '_space-individual_desc-rleftAmygdala_roi.nii']);
%roi_fns{2} = fullfile(anat_dir_preproc, ['sub-' sub '_space-individual_desc-rrightAmygdala_roi.nii']);
%roi_text  = 'bilateral amygdala';

[p, frm, rg, dim] = fmrwhy_util_readNifti(roi_fns{1});
roi_img = fmrwhy_util_createBinaryImg(p.nii.img, 0.1);
I_roi_distr = find(roi_img(:));

task_info.TR = 2;
task_info.onsets = [11; 31; 51; 71; 91; 111; 131; 151; 171; 191];
task_info.durations = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10];
task_info.precision = 1;

%fmrwhy_util_compareTSNR(tsnr_fns, mask_fn, roi_fns, roi_text)

%functional_fn = fullfile(func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srapreproc_bold.nii']);
functional_fn = fullfile(func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-srapreproc_bold.nii']);
fmrwhy_util_computeROImeasures(functional_fn, roi_img, task_info, tsnr_fn, 'stuff', 0, options)