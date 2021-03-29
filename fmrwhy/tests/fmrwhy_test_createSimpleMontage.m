
fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-qc/sub-001/func/sub-001_task-emotion_run-2_space-individual_tsnr.nii';
%fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/func/sub-001_task-motor_run-2_echo-2_desc-rapreproc_tsnr.nii';
[p, frm, rg, dim] = fmrwhy_util_readOrientNifti(fn);
img = p.nii.img;
overlaymontage = fmrwhy_util_createMontage(img, 9, 1, 'xxx', 'hot', 'on', 'max', []); colorbar;