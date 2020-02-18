function fmrwhy_batch_sliceTiming(functional_fn, saveAs_fn, defaults)

% Load required defaults

[d, f, e] = fileparts(functional_fn);
temp_functional_fn = fullfile(d, ['temp_' f e]);
copyfile(functional_fn, temp_functional_fn)
func_spm = spm_vol(temp_functional_fn);
Nt = numel(func_spm);
TR = defaults.TR;
N_slices = defaults.N_slices;

% Create cell array of scan names
scans = {};
for i = 1:Nt
    scans{i} = [temp_functional_fn ',' num2str(i)];
end

% Create SPM12 batch job
spm('defaults','fmri');
spm_jobman('initcfg');
slice_timing = struct;
slice_timing.matlabbatch{1}.spm.temporal.st.scans = {scans'};
slice_timing.matlabbatch{1}.spm.temporal.st.nslices = N_slices;
slice_timing.matlabbatch{1}.spm.temporal.st.tr = TR;
slice_timing.matlabbatch{1}.spm.temporal.st.ta = TR - TR/N_slices; % if 'to' is provided as array of slice times, this is ignored ==> set to zero
slice_timing.matlabbatch{1}.spm.temporal.st.so = 1:N_slices; % can also pass array of slice times in ms
slice_timing.matlabbatch{1}.spm.temporal.st.refslice = 1; % if 'to' is provided as array of slice times, refslice should be given in ms
slice_timing.matlabbatch{1}.spm.temporal.st.prefix = 'a';
% Run
spm_jobman('run',slice_timing.matlabbatch);
[d, f, e] = fileparts(temp_functional_fn);
atemp_functional_fn = fullfile(d, ['a' f e]);
[status, msg, msgID] = movefile(atemp_functional_fn, saveAs_fn);
delete(temp_functional_fn);