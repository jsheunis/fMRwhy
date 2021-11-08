function fmrwhy_batch_coregEst(anatomical_fn, template_fn, saveAs_fn)

    % First create a temporary copy of the anatomical_fn timeseries, since the
    % estimation process changes the header of the nifti files, and we want the original
    % image to remain unchanged
    [d, f, e] = fileparts(anatomical_fn);
    temp_anatomical_fn = fullfile(d, ['temp_' f e]);
    copyfile(anatomical_fn, temp_anatomical_fn);

    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    coreg_estimate = struct;
    % Ref
    coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.ref = {template_fn};
    % Source
    coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.source = {temp_anatomical_fn};
    % Other
    % coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.other = {};
    % Eoptions
    coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    % Run
    spm_jobman('run', coreg_estimate.matlabbatch);
    [status, msg, msgID] = movefile(temp_anatomical_fn, saveAs_fn);
