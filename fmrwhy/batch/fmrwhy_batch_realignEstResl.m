function fmrwhy_batch_realignEstResl(functional_fn, template_fn, saveAs_fn)

    [d, f, e] = fileparts(functional_fn);
    temp_functional_fn = fullfile(d, ['temp_' f e]);
    copyfile(functional_fn, temp_functional_fn);
    func_spm = spm_vol(temp_functional_fn);
    Nt = numel(func_spm);

    % Filenames for which to estimate 3D realignment parameters
    fns = {};
    if template_fn == 0
        for i = 1:Nt
            fns{i} = [temp_functional_fn ',' num2str(i)];
        end
    else
        fns{1} = [template_fn ',1'];
        for i = 1:Nt
            fns{i + 1} = [temp_functional_fn ',' num2str(i)];
        end
    end

    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    % Data
    realign_estimate_reslice = struct;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.data = {fns'};
    % Eoptions
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % register to first
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    % Roptions
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [1 0]; % images [2..n]
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    % Run
    spm_jobman('run', realign_estimate_reslice.matlabbatch);
    % Output
    output = struct;
    [d, f, e] = fileparts(temp_functional_fn);
    rtemp_functional_fn = fullfile(d, ['r' f e]);
    [status, msg, msgID] = movefile(rtemp_functional_fn, saveAs_fn);
    delete(temp_functional_fn);
