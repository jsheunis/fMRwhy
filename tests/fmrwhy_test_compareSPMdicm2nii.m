
options = fmrwhy_defaults;
% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';
% Loop through subjects, sessions, tasks, runs, etc
sub = '001';
% Loop through sessions, tasks, runs, etc
ses = '';
task = 'rest';
run = '1';
echo = '2';


% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% Update workflow params with subject functional derivative filenames
options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

%%
template_ts_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_echo-' options.template_echo '_bold.nii']);
template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);

%% 1 = spm, 2 nifti, 3 = dicm2nii

spm1 = spm_vol(template_ts_fn); 
img1 = spm_read_vols(spm1);
nii2 = nifti(template_ts_fn);
[p3, frm3, rg3, dim3] = fmrwhy_util_readNifti(template_ts_fn);
nii4 = nii_tool('load', template_ts_fn);


%% Plot same slice in images read in different ways
slice = 4;
figure;
subplot(141);
imagesc(rot90(squeeze(img1(:,:,slice,1)))); colorbar; title('spm read vols')
subplot(142);
imagesc(rot90(squeeze(nii2.dat(:,:,slice,1)))); colorbar; title('spm nifti')
subplot(143);
imagesc(rot90(squeeze(p3.nii.img(:,:,slice,1)))); colorbar; title('nii viewer - dicm2nii')
subplot(144);
imagesc(rot90(squeeze(nii4.img(:,:,slice,1)))); colorbar; title('nii tool - dicm2nii')

%%
% Plot
figure; subplot(121); imagesc(rot90(squeeze(p3.nii.img(:,:, 20,1)))); colorbar;
subplot(122); imagesc(rot90(squeeze(p4.nii.img(:,:,20)))); colorbar;



%%
fn1 = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_spmvol_image.nii']);
fn3 = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_dicm2n_image.nii']);

fn_converted = fullfile(options.sub_dir_preproc, 'func', 'converted.nii');
fn_converted1 = fullfile(options.sub_dir_preproc, 'func', 'converted_00001.nii');

%% save option 1: create duplicate nii structure, change hdr.dim field
% values to make it 3d, use nii_tool('save',)
new_nii = struct;
new_nii.hdr = nii4.hdr;
new_nii.img = nii4.img(:,:,:,1);
new_nii.hdr.dim(1) = 3; % set number of dimensions of image ==> 3D (not 4D)
new_nii.hdr.dim(5) = 1; % set the value of 4th dimension (time) to 1
new_nii.hdr.aux_file = '';
new_nii.hdr.file_name = fn3;
nii_tool('save', new_nii, fn3);

%% save option 2: convert 4d to 3d
nii_tool('save', nii4, fn_converted, true);

%% read in from saved 
%option 1
nii_saved1 = nii_tool('load', fn3);
% [psaved1, frm5, rg5, dim5] = fmrwhy_util_readNifti(fn3);

% option 2
nii_saved2 = nii_tool('load', fn_converted1);
% [psaved2, frm5, rg5, dim5] = fmrwhy_util_readNifti(fn3);

%% Plot same slice in images saved in different ways
slice = 30;
figure;
subplot(121);
imagesc(rot90(squeeze(nii_saved1.img(:,:,slice)))); colormap gray; colorbar; title('using original hdr')
subplot(122);
imagesc(rot90(squeeze(nii_saved2.img(:,:,slice)))); colormap gray; colorbar; title('4d to 3d conversion')


%%

%%
% Tsnr
[i,j,k,t]  = size(nii4.img);

img4d = single(nii4.img);
% img4dd = double(nii4.img);
img2d = reshape(img4d, i*j*k, t);
% img2dd = reshape(img4dd, i*j*k, t);
mn = mean(img2d, 2);
stddev = std(img2d, [], 2);
tsnr = mn./stddev;
tsnr_3d = reshape(tsnr, i, j, k);

% save using previous hdr
fn_tsnr = fullfile(options.sub_dir_preproc, 'func', 'tsnr3d.nii');
new_nii = struct;
new_nii.hdr = nii4.hdr;
new_nii.img = tsnr_3d;
new_nii.hdr.dim(1) = 3; % set number of dimensions of image ==> 3D (not 4D)
new_nii.hdr.dim(5) = 1; % set the value of 4th dimension (time) to 1
new_nii.hdr.aux_file = '';
new_nii.hdr.file_name = fn_tsnr;
nii_tool('save', new_nii, fn_tsnr);


%%



%%
%%
%%
%%
%%
%%

if ~exist(template_fn, 'file')
    disp(['Template funcional image does not exist yet. Creating now: ' template_fn]);
    functional0_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_echo-' options.template_echo '_bold.nii,1']);
    fmrwhy_util_saveNifti(template_fn, spm_read_vols(spm_vol(functional0_fn)), functional0_fn, 'Template functional volume', 0)
else
    disp(['Template funcional image exists: ' template_fn]);
end
options.template_fn = template_fn;