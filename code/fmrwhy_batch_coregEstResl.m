function fmrwhy_batch_coregEstResl(source_fn, reference_fn, other_fn, interp, saveAs_fn)

% First create a temporary copy of the source_fn timeseries, since the
% estimation process changes the header of the nifti files, and we want the original
% image to remain unchanged
[d, f, e] = fileparts(source_fn);
temp_source_fn = fullfile(d, ['temp_' f e]);
copyfile(source_fn, temp_source_fn)

spm('defaults','fmri');
spm_jobman('initcfg');
coreg_estwrite = struct;

% Ref
coreg_estwrite.matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {reference_fn};
% Source
coreg_estwrite.matlabbatch{1}.spm.spatial.coreg.estwrite.source = {temp_source_fn};
% Other
 coreg_estwrite.matlabbatch{1}.spm.spatial.coreg.estwrite.other = other_fn';
% Eoptions
coreg_estwrite.matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
coreg_estwrite.matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
coreg_estwrite.matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
coreg_estwrite.matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
% Roptions
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = interp;
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
% Run
spm_jobman('run',coreg_estwrite.matlabbatch);
%[status, msg, msgID] = movefile(temp_source_fn, saveAs_fn);
%delete(temp_source_fn);
