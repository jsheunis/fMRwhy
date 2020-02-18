function runSpmStats(jobs_dir, stats_dir, img_fn, Ndyn, Nregr, MR_fn)


%% specify design matrix in SPM.mat
design_stats_mat = [jobs_dir filesep 'model_specification.mat'];
design_stats = load(design_stats_mat);

design_stats.matlabbatch{1}.spm.stats.fmri_spec.dir = {stats_dir}; 
design_stats.matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
design_stats.matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
design_stats.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
design_stats.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

% if numel(design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess) > 1
%     design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(2:numel(design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess)) = [];
% end
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = {};
for i = 1:Ndyn
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans{i,1} = [img_fn ',' num2str(i) ];
end

design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.name = 'Word generation';
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.onset = [17;49;81;113;145;177];
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.duration = 16;
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {MR_fn};
% design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.tmod = 0;
% design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.pmod = {''};
% design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond.orth = 1;
% design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {};
% design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct;
% design_stats.matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
% design_stats.matlabbatch{1}.spm.stats.fmri_spec.bases.hrf = struct('derivs', [0 0]);
% design_stats.matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
% design_stats.matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
% design_stats.matlabbatch{1}.spm.stats.fmri_spec.volt = 0.8000;
% design_stats.matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
% design_stats.matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
% assume for now that rest of parameters stay the same as specified in
% initial job. this will change later or when requirements for script are
% different

cfg_util('run',design_stats.matlabbatch);
% spm_jobman('run',design_stats);

%% Estimate model
% either use spm_spm() or saved spm job method (TODO: test both)
% SPM = load([stats_dir filesep 'SPM.mat']); 
% unsure about this methos, maybe it is only used for second level analysis
% spm_spm(SPM);

model_estimation_mat = [jobs_dir filesep 'model_estimation.mat'];
model_estimation = load(model_estimation_mat);
model_estimation.matlabbatch{1}.spm.stats.fmri_est.spmmat = {[stats_dir filesep 'SPM.mat']};
model_estimation.matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
cfg_util('run',model_estimation.matlabbatch);


%% Create and run contrasts
contrast_mat = [jobs_dir filesep 'contrast_setup.mat'];
contrast = load(contrast_mat);
contrast.matlabbatch{1}.spm.stats.con.spmmat = {[stats_dir filesep 'SPM.mat']};
contrast.matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Task';
weights = zeros(1, 1+Nregr+1); weights(1) = 1;
contrast.matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = weights;
contrast.matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
contrast.matlabbatch{1}.spm.stats.con.delete = 0;
cfg_util('run',contrast.matlabbatch);

%% Run results job
spm('defaults','FMRI');
stats_results_mat = [jobs_dir filesep 'stats_results.mat'];
stats_results = load(stats_results_mat);
stats_results.matlabbatch{1}.spm.stats.results.spmmat = {[stats_dir filesep 'SPM.mat']};
% assume for now that rest of parameters stay the same as specified in
% initial job. this will change later or when requirements for script are
% different
cfg_util('run',stats_results.matlabbatch);







