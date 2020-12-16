function fmrwhy_util_convert3Dto4D(fns, TR, out_fn)

% Initialize
spm('defaults','fmri');
spm_jobman('initcfg');
% Setup structure
convert3dto4d = struct;
convert3dto4d.matlabbatch{1}.spm.util.cat.vols = fns';
convert3dto4d.matlabbatch{1}.spm.util.cat.name = out_fn;
convert3dto4d.matlabbatch{1}.spm.util.cat.dtype = 0; % same as input dtypes
convert3dto4d.matlabbatch{1}.spm.util.cat.RT = 2;
% Run
spm_jobman('run',convert3dto4d.matlabbatch);