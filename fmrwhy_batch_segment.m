function fmrwhy_batch_segment(anatomical_fn, spm_dir, saveAs_fns, saveAsTransforms_fns)

spm('defaults','fmri');
spm_jobman('initcfg');
segmentation = struct;
% Channel
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.vols = {anatomical_fn};
% Tissue
for t = 1:6
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm = {[spm_dir filesep 'tpm' filesep 'TPM.nii,' num2str(t)]};
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus = t-1;
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [1 0];
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).warped = [0 0];
end
segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
% Warp
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.write=[1 1];
% Run
spm_jobman('run',segmentation.matlabbatch);

% Save filenames
[d, fn, ext] = fileparts(anatomical_fn);
% Save segmentation filenames: 1=gm, 2=wm, 3=csf, 4=bone, 5=soft_tissue, 6=air
for i = 1:6
    cfn = fullfile(d, ['c' num2str(i) fn ext]);
    [status, msg, msgID] = movefile(cfn, saveAs_fns{i});
end
% Save transform filenames
forward_transformation = fullfile(d, ['y_' fn ext]);
inverse_transformation = fullfile(d, ['iy_' fn ext]);
[status, msg, msgID] = movefile(forward_transformation, saveAsTransforms_fns{1});
[status, msg, msgID] = movefile(inverse_transformation, saveAsTransforms_fns{2});