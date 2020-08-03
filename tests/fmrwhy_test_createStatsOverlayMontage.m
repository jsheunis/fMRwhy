% Template
%template_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-coregEstResl_T1w.nii';
template_tsnr_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-multiecho/sub-001/func/sub-001_task-rest_run-1_echo-2_desc-rapreproc_tsnr.nii';
template_perc_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-multiecho/sub-001/func/sub-001_task-rest_run-1_echo-2_desc-rapreproc_tsnr.nii';

[p, frm, rg, dim] = fmrwhy_util_readOrientNifti(template_tsnr_fn);
template_img = p.nii.img;

% Stats

%tmap = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-motor_run-1/spmT_0001.nii';
%tmap_clusters = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-motor_run-1/spmT_0001_nary_clusters.nii';
tmap = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-emotion_run-1/spmT_0001.nii';
tmap_clusters = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-emotion_run-1/spmT_0001_nary_clusters.nii';
[ptmap, frm3, rg3, dim3] = fmrwhy_util_readOrientNifti(tmap);
[ptmapc, frm3, rg3, dim3] = fmrwhy_util_readOrientNifti(tmap_clusters);
stats_img = fmrwhy_util_maskImage(double(ptmap.nii.img), double(ptmapc.nii.img));
stats_img = [];

% ROIs
roi_lmotor_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-rleftMotor_roi.nii';
roi_bamygdala_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-rBilateralAmygdala_roi.nii';
roi_fns = {roi_lmotor_fn, roi_bamygdala_fn};
[p1, frm1, rg1, dim1] = fmrwhy_util_readOrientNifti(roi_lmotor_fn);
[p2, frm2, rg2, dim2] = fmrwhy_util_readOrientNifti(roi_bamygdala_fn);
%roi_img = {};
%roi_img{1} = fmrwhy_util_createBinaryImg(p1.nii.img, 0.1);
%roi_img{2} = fmrwhy_util_createBinaryImg(p2.nii.img, 0.1);
%roi_img = fmrwhy_util_createBinaryImg(p1.nii.img, 0.1);
roi_img = fmrwhy_util_createBinaryImg(p2.nii.img, 0.1);

roi_img = {};
I_roi = {};
overlay_img = zeros(dim);
for i = 1:numel(roi_fns)
    [p, ~, ~, ~] = fmrwhy_util_readOrientNifti(roi_fns{i});
    roi_img{i} = fmrwhy_util_createBinaryImg(p.nii.img, 0.1);
    I_roi{i} = find(roi_img{i}(:));
    overlay_img = overlay_img | roi_img{i};
end





% Parameters
saveAs_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-multiecho/sub-001/blabla.png';
saveAs_fn = [];
columns = 9;
rotate = 1;
str = '';
%clrmp = 'hot';
clrmp = 'parula';
visibility = 'on';
shape = 'max';
%cxs = [0 250];
cxs = [0 300];
stats_clrmp = [];
clrbar = true;
%roi_rgbcolors = [148, 239, 255];
roi_rgbcolors = [148, 239, 255];

% Call function
%output = fmrwhy_util_createStatsOverlayMontage(template_img, stats_img, roi_img, columns, rotate, str, clrmp, visibility, shape, cxs, stats_clrmp, roi_rgbcolors, saveAs_fn)
output = fmrwhy_util_createStatsOverlayMontage(template_img, stats_img, overlay_img, columns, rotate, str, clrmp, visibility, shape, cxs, stats_clrmp, roi_rgbcolors, clrbar, saveAs_fn)