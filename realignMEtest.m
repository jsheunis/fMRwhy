
% Load required defaults
preproc_dir = '/Users/jheunis/Desktop/realignMEtest';
sub = 'sub-pilot';
template_task = 'rest';
template_run = 1;
template_echo = 2;

functional_fn = fullfile(preproc_dir, sub, 'func', [sub '_task-' template_task '_run-' num2str(template_run) '_echo-' num2str(template_echo) '.nii']);
functional0_fn = [functional_fn ',1'];
template_vol = fullfile(preproc_dir, sub, 'func', [sub '_task-' template_task '_run-' num2str(template_run) '_echo-' num2str(template_echo) '_template.nii']);
rtme_util_saveNifti(functional0_fn, spm_read_vols(spm_vol(functional0_fn)), template_vol, 'Template functional volume')

%%

task = 'rest';
run = 1;
reference_echo = 2;
N_vol = 210;

% Create cell array of scan names to realign template echo timeseries
functional_fn = fullfile(preproc_dir, sub, 'func', [sub '_task-' task '_run-' num2str(run) '_echo-' num2str(reference_echo) '.nii']);
scans = {};
scans{1} = template_vol;
for i = 2:N_vol+1
    scans{i} = [functional_fn ',' num2str(i-1)];
end

% Realign (estimate and reslice) template echo timeseries
spm('defaults','fmri');
spm_jobman('initcfg');
realign_estimate_reslice = struct;
% Data
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.data={scans'};
% Eoptions
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % register to first
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
% Roptions
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [1 0]; % images [2..n]
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
% Run
spm_jobman('run',realign_estimate_reslice.matlabbatch);
output.rfunctional_fn = [d filesep 'r' fn ext];

%%
% Access realignment (head motion) parameters
[d, fn, ext] = fileparts(template_vol);
output.HMP_temp_fn = [d filesep 'rp_' fn '.txt'];
HMP = load(output.HMP_temp_fn);

HMP(1,:) = [];
output.HMP_fn = fullfile(d, ['rp_' sub '_task-' task '_run-' num2str(run) '.txt']);
dlmwrite(output.HMP_fn, HMP, 'delimiter', '\t', 'precision', '%1.7e')

delete(output.HMP_temp_fn)

%%

e = 1;
functional_fn = fullfile(preproc_dir, sub, 'func', [sub '_task-' task '_run-' num2str(run) '_echo-' num2str(e) '.nii']);
rfunctional_fn = fullfile(preproc_dir, sub, 'func', ['r' sub '_task-' task '_run-' num2str(run) '_echo-' num2str(e) '.nii']);
functional_spm = spm_vol(functional_fn);
functional_img = spm_read_vols(functional_spm);
rfunctional_spm = functional_spm;
for i = 1:N_vol
    currentVol = functional_spm(i);
    currentImg = functional_img(:,:,:,i);
    Pm = zeros(12,1);
    Pm(1:6) = HMP(i, :);
    orig_mat = currentVol.mat;
    rigid_mat = spm_matrix(Pm, 'T*R');
    trans_mat = rigid_mat * orig_mat;
    rfunctional_spm(i).mat = trans_mat;
    rfunctional_spm(i).fname = rfunctional_fn;
    rfunctional_spm(i).private.dat.fname = rfunctional_fn;
    spm_write_vol(rfunctional_spm(i),functional_img(:,:,:,i));
end
disp('fnishe')
%%
scans = {};
scans{1} = template_vol;
for i = 2:N_vol+1
    scans{i} = [rfunctional_fn ',' num2str(i-1)];
end

% Reslice echo timeseries
spm('defaults','fmri');
spm_jobman('initcfg');
realign_reslice = struct;
% Data
realign_reslice.matlabbatch{1}.spm.spatial.realign.write.data=scans';
% Roptions
realign_reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.which = [1 0]; % images [2..n]
realign_reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
realign_reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
realign_reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
realign_reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
% Run
spm_jobman('run',realign_reslice.matlabbatch);
%%
[d, fn, ext] = fileparts(rfunctional_fn);
rrfunctional_fn = [d filesep 'r' fn ext];

%%
delete(rfunctional_fn)
[status, msg, msgID] = movefile(rrfunctional_fn, rfunctional_fn);
if status == 0
    disp(msg)
end
