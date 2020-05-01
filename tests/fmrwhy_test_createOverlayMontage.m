

template_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/func/sub-001_task-rest_run-1_space-individual_bold.nii';
anat_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_T1w.nii';
ranat_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-coregEstResl_T1w.nii';

overlay_fn = {};
overlay_fn{1} = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-rLeftAmygdala_roi.nii';
overlay_fn{2} = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-rrightAmygdala_roi.nii';
overlay_fn{3} = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-rLeftMotor_roi.nii';
overlay_fn{4} = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-rrightMotor_roi.nii';

%[p1, frm1, rg1, dim1] = fmrwhy_util_readNifti(template_fn);
%[p2, frm2, rg2, dim2] = fmrwhy_util_readNifti(ranat_fn);
%
%template_img = p2.nii.img;
%overlay_img = {};
%for i=1:numel(overlay_fn)
%    [p3, frm3, rg3, dim3] = fmrwhy_util_readNifti(overlay_fn{i});
%    overlay_img{i} = fmrwhy_util_createBinaryImg(p3.nii.img, 0.1);
%end
%
%rotate = 1;
%columns = 9;
%rotate = 1;
%str = '';
%clrmp = 'gray';
%visibility = 'on';
%shape = 'max';
%saveAs_fn = 'overlaytest';
%
%saveAs_fn = {};
%saveAs_fn{1} = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-qc/sub-001/anat/sub-001_space-individual_desc-leftAmygdala_roi_montage.png';
%saveAs_fn{2} = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-qc/sub-001/anat/sub-001_space-individual_desc-rightAmygdala_roi_montage.png';
%saveAs_fn{3} = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-qc/sub-001/anat/sub-001_space-individual_desc-leftMotor_roi_montage.png';
%saveAs_fn{4} = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-qc/sub-001/anat/sub-001_space-individual_desc-rightMotor_roi_montage.png';
%
%for i=1:4
%    output = fmrwhy_util_createOverlayMontage(template_img, overlay_img{i}, columns, rotate, str, clrmp, visibility, shape, saveAs_fn{i});
%end
%
%
%%montage1 = fmrwhy_util_createMontage(p1.nii.img, 9, rotate, 'test', 'gray', 'off', 'max');
%%montage2 = fmrwhy_util_createMontage(p2.nii.img(:,:,181:214), 9, rotate, 'test', 'gray', 'off', 'max');

koek = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-WM_mask.nii';
[p3, frm3, rg3, dim3] = fmrwhy_util_readOrientNifti(koek);
roi_img = fmrwhy_util_createBinaryImg(p3.nii.img, 0.1);
fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/func/sub-001_task-rest_run-1_echo-2_desc-rapreproc_tsnr.nii';
template_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/func/sub-001_task-rest_run-1_space-individual_bold.nii';
[p, frm, rg, dim] = fmrwhy_util_readOrientNifti(template_fn);
%montage_template = fmrwhy_util_createMontage(p.nii.img, 9, 1, 'Template volume', 'hot', 'on', 'max', [0 100]); colorbar;
overlaymontage = fmrwhy_util_createOverlayMontage(p.nii.img, roi_img, 9, 1, 'xxx', 'gray', 'on', 'max', [], 0);