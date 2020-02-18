function fmrwhy_batch_smooth(functional_fn, saveAs_fn, fwhm)

[d, f, e] = fileparts(functional_fn);
temp_functional_fn = fullfile(d, ['temp_' f e]);
copyfile(functional_fn, temp_functional_fn)
func_spm = spm_vol(temp_functional_fn);
Nt = numel(func_spm);

% Create cell array of scan names
scans = {};
for i = 1:Nt
    scans{i} = [temp_functional_fn ',' num2str(i)];
end

% Create SPM12 batch job
spm('defaults','fmri');
spm_jobman('initcfg');
smooth = struct;
smooth.matlabbatch{1}.spm.spatial.smooth.data = scans';
smooth.matlabbatch{1}.spm.spatial.smooth.fwhm = fwhm;
smooth.matlabbatch{1}.spm.spatial.smooth.dtype = 0;
smooth.matlabbatch{1}.spm.spatial.smooth.im = 0;
smooth.matlabbatch{1}.spm.spatial.smooth.prefix = 's';
% Run
spm_jobman('run',smooth.matlabbatch);

[d, f, e] = fileparts(temp_functional_fn);
stemp_functional_fn = fullfile(d, ['s' f e]);
[status, msg, msgID] = movefile(stemp_functional_fn, saveAs_fn);
delete(temp_functional_fn);