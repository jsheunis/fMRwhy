 function viewSpmStatsResults(jobs_dir, stats_dir)
spm('defaults','FMRI');
stats_results_mat = [jobs_dir filesep 'stats_results.mat'];
stats_results = load(stats_results_mat);
stats_results.matlabbatch{1}.spm.stats.results.spmmat = {[stats_dir filesep 'SPM.mat']};
% assume for now that rest of parameters stay the same as specified in
% initial job. this will change later or when requirements for script are
% different
cfg_util('run',stats_results.matlabbatch);