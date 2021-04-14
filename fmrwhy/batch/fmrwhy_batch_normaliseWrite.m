function fmrwhy_batch_normaliseWrite(toTransform_fns, transformation_fn, template_fn, saveAs_fns)

    % toTransform_fns = cell array of filenames; TODO, check orientation (column array or row array)

    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    normalize_write = struct;
    % Ref
    normalize_write.matlabbatch{1}.spm.spatial.normalise.write.subj.def = {transformation_fn};
    % Data
    normalize_write.matlabbatch{1}.spm.spatial.normalise.write.subj.resample = toTransform_fns';
    % Write options
    normalize_write.matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78, -112, -70; 78, 76, 85];
    normalize_write.matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2, 2, 2];
    normalize_write.matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    normalize_write.matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    % Run
    spm_jobman('run', normalize_write.matlabbatch);

    % Save resliced filenames
    for i = 1:numel(toTransform_fns)
        [d, fn, ext] = fileparts(toTransform_fns{i});
        wfn = fullfile(d, ['w' fn ext]);
        [status, msg, msgID] = movefile(wfn, saveAs_fns{i});
    end
