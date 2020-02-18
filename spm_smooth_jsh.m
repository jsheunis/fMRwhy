function output = spm_smooth_jsh(functional_fn, fwhm)

func_spm = spm_vol(functional_fn);
Nt = numel(func_spm);
smooth = struct;
% Data
fns={};
for i = 1:Nt
    fns{i} = [functional_fn ',' num2str(i) ];
end
smooth.matlabbatch{1}.spm.spatial.smooth.data = fns';
% Other
smooth.matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
smooth.matlabbatch{1}.spm.spatial.smooth.dtype = 0;
smooth.matlabbatch{1}.spm.spatial.smooth.im = 0;
smooth.matlabbatch{1}.spm.spatial.smooth.prefix = 's';
% Run
cfg_util('run',smooth.matlabbatch);
[d, f, e] = fileparts(functional_fn);
output.sfunctional_fn = [d filesep 's' f e];