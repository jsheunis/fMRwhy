function fmrwhy_batch_specify2ndlevel(stats_dir, con_fns, params)

spm('defaults','fmri');
spm_jobman('initcfg');

% SETUP BATCH JOB STRUCTURE
design_stats = struct;
% dir
design_stats.matlabbatch{1}.spm.stats.factorial_design.dir = {stats_dir};
% des.t1.scans
design_stats.matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = con_fns';
% cov
design_stats.matlabbatch{1}.spm.stats.factorial_design.cov = struct([]);
% multi_cov
design_stats.matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct([]);
% masking
design_stats.matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
design_stats.matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
design_stats.matlabbatch{1}.spm.stats.factorial_design.masking.em = {[]};
% globalc
design_stats.matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
% globalm
design_stats.matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
design_stats.matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% RUN BATCH JOB
spm_jobman('run',design_stats.matlabbatch);