function fmrwhy_batch_realignResl(reslice_scans)

    %% ---------
    %% ---------
    %% ---------
    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    reslice = struct;
    % Ref
    reslice.matlabbatch{1}.spm.spatial.realign.write.data = reslice_scans';
    % Roptions
    reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.which = [1 0]; % images [2..n]
    reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
    reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
    reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 0;
    reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'rr';
    % Run
    spm_jobman('run', reslice.matlabbatch);
