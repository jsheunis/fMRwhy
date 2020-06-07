% Template
template_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-coregEstResl_T1w.nii';
[p, frm, rg, dim] = fmrwhy_util_readOrientNifti(template_fn);
template_img = p.nii.img;

% Stats

%tmap = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-motor_run-1/spmT_0001.nii';
%tmap_clusters = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-motor_run-1/spmT_0001_nary_clusters.nii';
tmap = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-emotion_run-1/spmT_0001.nii';
tmap_clusters = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-emotion_run-1/spmT_0001_nary_clusters.nii';
[ptmap, frm3, rg3, dim3] = fmrwhy_util_readOrientNifti(tmap);
[ptmapc, frm3, rg3, dim3] = fmrwhy_util_readOrientNifti(tmap_clusters);
stats_img = fmrwhy_util_maskImage(double(ptmap.nii.img), double(ptmapc.nii.img));

% ROIs
roi_lmotor_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-rleftMotor_roi.nii';
roi_bamygdala_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-rBilateralAmygdala_roi.nii';
[p1, frm1, rg1, dim1] = fmrwhy_util_readOrientNifti(roi_lmotor_fn);
[p2, frm2, rg2, dim2] = fmrwhy_util_readOrientNifti(roi_bamygdala_fn);
%roi_img = {};
%roi_img{1} = fmrwhy_util_createBinaryImg(p1.nii.img, 0.1);
%roi_img{2} = fmrwhy_util_createBinaryImg(p2.nii.img, 0.1);
%roi_img = fmrwhy_util_createBinaryImg(p1.nii.img, 0.1);
roi_img = fmrwhy_util_createBinaryImg(p2.nii.img, 0.1);

% Parameters
saveAs_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-motor_run-1/blabla.png';
columns = 9;
rotate = 1;
str = '';
clrmp = 'gray';
visibility = 'off';
shape = 'max';
cxs = [];
stats_clrmp = 'hot';
roi_rgbcolors = [148, 239, 255];

% Call function
output = fmrwhy_util_createStatsOverlayMontage(template_img, stats_img, roi_img, columns, rotate, str, clrmp, visibility, shape, cxs, stats_clrmp, roi_rgbcolors, saveAs_fn)