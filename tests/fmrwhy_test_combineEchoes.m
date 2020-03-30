% Load/create required parameters
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, []);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Loop through subjects, sessions, tasks, runs, etc
sub = '001';

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% Loop through sessions, tasks, runs, etc
ses = '';
task = 'motor';
run = '1';
echo = '2';

%% Update workflow params with subject functional derivative filenames
options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

% Grab/construct parameters for multi-echo combination

disp('Concatenating functional data')
TE = options.TE;
[p, frm, rg, template_dim] = fmrwhy_util_readNifti(options.template_fn);
Nt = options.Nscans;
sz = [template_dim Nt numel(TE)];
func_data = zeros(sz);
for e = 1:numel(TE)
    echo = num2str(e);
    rafunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_bold.nii']);
    [p, frm, rg, dim] = fmrwhy_util_readNifti(rafunctional_fn);
    func_data(:,:,:,:,e) = p.nii.img;
end

%disp('Loading mask')
%masks = fmrwhy_util_loadMasks(bids_dir, sub);
%mask = masks.brain_mask_3D;

disp('Loading weight images')
t2star_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_desc-MEparams_t2star.nii']);
[p, frm, rg, dim] = fmrwhy_util_readNifti(t2star_fn);
t2star_img = p.nii.img;

tsnr_data = zeros([template_dim numel(TE)]);
for e = 1:numel(TE)
    echo = num2str(e);
    tsnr_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_echo-' echo '_desc-rapreproc_tsnr.nii']);
    [p, frm, rg, dim] = fmrwhy_util_readNifti(tsnr_fn);
    tsnr_data(:,:,:,e) = p.nii.img;
end

% data      - 5D [4D x E] timeseries or 4D [3D x E] volume data to be combined_data
% TE        - vector of echo times (ms) in order of acquisition, e.g. [14 28 42]

TE = options.TE;
data = squeeze(func_data(:,:,:,10,:)); %echo 1 timepoint 10
tic;
combined_data_t2s = fmrwhy_me_combineEchoes(data, TE, 0, 1, t2star_img);
toc;
%combined_dataAll_t2s = fmrwhy_me_combineEchoes(func_data, TE, 0, 1, t2star_img);
%combined_data_tsnr = fmrwhy_me_combineEchoes(data, TE, 0, 2, tsnr_data);
%combined_dataAll_tsnr = fmrwhy_me_combineEchoes(func_data, TE, 0, 2, tsnr_data);
%combined_dataAll_TE = fmrwhy_me_combineEchoes(func_data, TE, 0, 3, TE);
%combined_data_TE = fmrwhy_me_combineEchoes(data, TE, 0, 3, TE);

img1 = squeeze(func_data(:,:,:,10,1)); %echo 1 timepoint 10
img2 = squeeze(func_data(:,:,:,10,2)); %echo 2 timepoint 10
img3 = squeeze(func_data(:,:,:,10,3)); %echo 3 timepoint 10
%img4 = squeeze(combined_dataAll_t2s(:,:,:,10)); % combined data timepoint 10
img4 = combined_data_t2s; % combined data timepoint 10
%img4 = squeeze(combined_dataAll_tsnr(:,:,:,10)); % combined data timepoint 10
%img4 = squeeze(combined_dataAll_TE(:,:,:,10)); % combined data timepoint 10
columns = 9;
rotate = 1;
str = '';
clrmp = 'gray';
visibility = 'on';
shape = 'max';
montage1 = fmrwhy_util_createMontage(img1, columns, rotate, str, clrmp, visibility, shape)
montage2 = fmrwhy_util_createMontage(img2, columns, rotate, str, clrmp, visibility, shape)
montage3 = fmrwhy_util_createMontage(img3, columns, rotate, str, clrmp, visibility, shape)
montage4 = fmrwhy_util_createMontage(img4, columns, rotate, str, clrmp, visibility, shape)

%combined_dataAll_tsnr = fmrwhy_me_combineEchoes(func_data, TE, mask, method, weight_data);
%combined_dataAll_t2s = fmrwhy_me_combineEchoes(func_data, TE, mask, method, weight_data);
%combined_dataAll_tsnr = fmrwhy_me_combineEchoes(func_data, TE, mask, method, weight_data);
%montage = fmrwhy_util_createMontage(img, columns, rotate, str, clrmp, visibility, shape)
