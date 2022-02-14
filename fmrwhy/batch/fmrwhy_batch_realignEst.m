function output = fmrwhy_batch_realignEst(functional_fn, template_fn)

    % First create a temporary copy of the functional timeseries, since the
    % estimation process changes the header of the nifti files, and we want the original
    % timeseries to remain unchanged
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
    realign_estimate = struct;
    realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.data = {fns'};
    % Eoptions
    realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
    realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
    realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
    realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 0; % register to first
    realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
    realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
    realign_estimate.matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';
    % Run
    spm_jobman('run', realign_estimate.matlabbatch);
    % Output
    output = struct;
    [d, f, e] = fileparts(temp_functional_fn);
    if template_fn == 0
        output.mp_fn = fullfile(d, ['rp_' f '.txt']);
    else
        [dref, fref, eref] = fileparts(template_fn);
        output.mp_fn = [dref filesep 'rp_' fref '.txt'];
    end
    output.MP = load(output.mp_fn);
    [rows, cols] = size(output.MP);
    if rows == Nt + 1
        output.MP(1, :) = [];
        dlmwrite(output.mp_fn, output.MP, 'delimiter', '\t', 'precision', '%1.7e');
    end

    delete(temp_functional_fn);
