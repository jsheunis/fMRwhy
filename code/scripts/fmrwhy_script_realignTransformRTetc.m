
sub_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/';
template_fn = fullfile(sub_dir, 'sub-010_task-rest_run-1_space-individual_bold.nii');
orig_bold_fn = fullfile(sub_dir, 'sub-010_task-emotion_run-1_echo-2_bold.nii');

%rspm_fn = fullfile(sub_dir, 'rsub-010_task-emotion_run-1_echo-2_bold.nii');
%rfmrwhy_fn = fullfile(sub_dir, 'sub-010_task-emotion_run-1_echo-2_desc-rpreproc_bold.nii');
%rspm_rfmrwhy_fn = fullfile(sub_dir, 'rsub-010_task-emotion_run-1_echo-2_desc-rpreproc_bold.nii');


% Step 1 - Estimate realignment parameters
functional_fn = orig_bold_fn;
%realign_measures = fmrwhy_batch_realignEst(functional_fn, template_fn);
%
%% Step 2 - Estimate and reslice
%functional_fn = orig_bold_fn;
%rfmrwhy_estResl_fn = fullfile(sub_dir, 'rfmrwhy_estResl.nii');
%saveAs_fn = rfmrwhy_estResl_fn;
%fmrwhy_batch_realignEstResl(functional_fn, template_fn, saveAs_fn)

% Step 3 - Apply MPs
toTransform_fn = orig_bold_fn;
motion_params = load(fullfile(sub_dir, 'rp_sub-010_task-rest_run-1_space-individual_bold_step1.txt'));
rfmrwhy_applyTransform_fn = fullfile(sub_dir, 'rt2fmrwhy_applyTransform.nii');
saveAs_fn = rfmrwhy_applyTransform_fn;
[tfunctional_spm, reslVol] = fmrwhy_util_applyTransformRT2(toTransform_fn, motion_params, template_fn, saveAs_fn, 2);

%
%functional_spm = spm_vol(orig_bold_fn);
%functional_img = spm_read_vols(functional_spm);
%N_vol = numel(functional_spm);
%origVol = {};
%% For each volume
%for i = 1:N_vol
%    origVol{i} = functional_img(:,:,:,i);
%end
%
%new_spm = spm_vol(orig_bold_fn);
%rtr2fmrwhy_applyTransform_fn = fullfile(sub_dir, 'rtr2fmrwhy_applyTransform.nii');
%
%for i = 1:N_vol
%    rfunctional_img(:,:,:, i) = reslVol{i};
%    new_spm(i).fname = rtrfmrwhy_applyTransform_fn;
%    new_spm(i).private.dat.fname = rtrfmrwhy_applyTransform_fn;
%    spm_write_vol(new_spm(i),rfunctional_img(:,:,:,i));
%end
%
%







% Movement parameters of SPM GUI realignment of original bold timeseries
%MPbold_fn = fullfile(sub_dir, 'MPbold.tsv');
%MPbold_struct = tdfread(MPbold_fn);
%MPbold_mat = struct2array(MPbold_struct);
%figure; plot(MPbold_mat);

% Movement parameters of SPM GUI realignment of fMRwhy realigned timeseries
%MPrapreproc_fn = fullfile(sub_dir, 'MPrapreproc.tsv');
%MPrapreproc_struct = tdfread(MPrapreproc_fn);
%MPrapreproc_mat = struct2array(MPrapreproc_struct);
%figure; plot(MPrapreproc_mat);

% Step 1: estimate realignment parameters from original bold timeseries with fMRwhy
%realign_measures = fmrwhy_batch_realignEst(orig_bold_fn, template_fn);
%MPbold_fmrwhy_txt = fullfile(sub_dir, 'MPbold_fmrwhy.txt');
%col_names = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
%dlmwrite(MPbold_fmrwhy_txt, realign_measures.MP, 'delimiter', '\t', 'precision', '%1.7e')
%MPbold_fmrwhy_tsv = fmrwhy_util_saveAsTSV(MPbold_fmrwhy_txt, col_names);
%MPbold_fmrwhy_fn = fullfile(sub_dir, 'MPbold_fmrwhy.tsv');
%MPbold_fmrwhy_struct = tdfread(MPbold_fmrwhy_fn);
%MPbold_fmrwhy_mat = struct2array(MPbold_fmrwhy_struct);
%figure; plot(MPrapreproc_mat);

% Result: this yields MPs that are the same as MPrapreproc ???

% Step 2: apply MP















