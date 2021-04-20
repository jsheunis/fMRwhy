function fmrwhy_batch_review1stlevel(stats_dir, params)

    % SETUP BATCH JOB STRUCTURE
    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    review = struct;
    review.matlabbatch{1}.spm.stats.review.spmmat = {[stats_dir filesep 'SPM.mat']};
    review.matlabbatch{1}.spm.stats.review.display.(params.display) = 1; % 'matrix', 'covariance', 'orth'
    review.matlabbatch{1}.spm.stats.review.print = params.print; % 'jpg', 'png',

    % RUN BATCH JOB
    spm_jobman('run', review.matlabbatch);
