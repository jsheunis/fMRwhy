function fmrwhy_batch_normaliseEstWrite(reference_fn, toTransform_fns, tpm_fn, saveAs_fns)

    % toTransform_fns = cell array of filenames; could be 3D or 4D niftis
    resample_fns = {};
    j = 0;
    for i = 1:numel(toTransform_fns)

        [d, f, e] = fileparts(toTransform_fns{i});
        temp_fn{i} = fullfile(d, ['temp_' f e]);
        copyfile(toTransform_fns{i}, temp_fn{i});
        func_spm = spm_vol(temp_fn{i});
        Nt = numel(func_spm);
        if Nt > 1
            for n = 1:Nt
                j = j+1;
                resample_fns{j} = [temp_fn{i} ',' num2str(n)];
            end
        else
            j = j+1;
            resample_fns{j} = [temp_fn{i} ',1'];
        end
    end

    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    normalize_estwrite = struct;
    % Ref
    normalize_estwrite.matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {reference_fn};
    % Data
    normalize_estwrite.matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = resample_fns';
    % Estimate options
    normalize_estwrite.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    normalize_estwrite.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    normalize_estwrite.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {tpm_fn};
    normalize_estwrite.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    normalize_estwrite.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    normalize_estwrite.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    normalize_estwrite.matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    % Write options
    normalize_estwrite.matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78, -112, -70; 78, 76, 85];
    normalize_estwrite.matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2, 2, 2];
    normalize_estwrite.matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    normalize_estwrite.matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
    % Run
    spm_jobman('run', normalize_estwrite.matlabbatch);

    % Save resliced filenames
    for i = 1:numel(toTransform_fns)
        [d, fn, ext] = fileparts(toTransform_fns{i});
        wfn = fullfile(d, ['w' fn ext]);
        [status, msg, msgID] = movefile(wfn, saveAs_fns{i});
    end