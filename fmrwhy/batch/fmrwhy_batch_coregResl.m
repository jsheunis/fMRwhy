function fmrwhy_batch_coregResl(reslice_fns, template_fn, saveAs_fns)

    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    reslice = struct;
    % Ref
    reslice.matlabbatch{1}.spm.spatial.coreg.write.ref = {template_fn};
    % Source
    reslice.matlabbatch{1}.spm.spatial.coreg.write.source = reslice_fns';
    % Roptions
    reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    % Run
    spm_jobman('run', reslice.matlabbatch);

    % Save resliced filenames
    for i = 1:numel(reslice_fns)
        [d, fn, ext] = fileparts(reslice_fns{i});
        rfn = fullfile(d, ['r' fn ext]);
        [status, msg, msgID] = movefile(rfn, saveAs_fns{i});
    end
