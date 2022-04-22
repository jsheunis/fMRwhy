function fmrwhy_batch_specify1stlevel(stats_dir, functional_fn, multi_reg_fn, sess_params, timing_params)

    spm('defaults', 'fmri');
    spm_jobman('initcfg');

    if iscell(functional_fn)
        if numel(functional_fn) ~= numel(multi_reg_fn) || numel(functional_fn) ~= numel(sess_params)
            error('For multirun specifications, arguments functional_fn, multi_reg_fn, and sess_params need to have the same number of elements')
        end
    else
        functional_fn = {functional_fn};
        multi_reg_fn = {multi_reg_fn};
        sess_params = {sess_params};
    end
    
    % SETUP BATCH JOB STRUCTURE
    design_stats = struct;
    % dir
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.dir = {stats_dir};
    % timing
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.timing.units = timing_params.units;
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.timing.RT = timing_params.RT;
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = timing_params.fmri_t;
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = timing_params.fmri_t0;
    % sess
    for s = 1:numel(functional_fn)
        func_fn = functional_fn{s};
        func4D_spm = spm_vol(func_fn);
        Nt = numel(func4D_spm);
        design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans = {};
        for i = 1:Nt
            design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(s).scans{i, 1} = [func_fn ',' num2str(i)];
        end
        design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(s).cond = sess_params{s}.cond;
        design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(s).multi = {''};
        design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(s).regress = {''};
        design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(s).multi_reg = {multi_reg_fn{s}};
        design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess(s).hpf = 128;
    end
    % fact
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.fact = {''};
    % bases
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.bases.hrf = struct('derivs', [0 0]);
    % volt
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    % global
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    % mthresh
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8000;
    % mask
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    % cvi
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

    % RUN BATCH JOB
    spm_jobman('run', design_stats.matlabbatch);
