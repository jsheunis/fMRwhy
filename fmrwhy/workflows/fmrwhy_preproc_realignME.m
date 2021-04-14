function output = rtme_preproc_realignME(sub, task, run, reference_echo, prefix, previous_steps_prefix, defaults)

    % Load required defaults
    TR = defaults.TR;
    N_slices = defaults.N_slices;
    N_echoes = defaults.N_echoes;
    N_vol = defaults.N_vol;
    preproc_dir = defaults.preproc_dir;
    template_vol = fullfile(preproc_dir, sub, 'func', [sub '_task-' template_task '_run-' num2str(template_run) '_echo-' num2str(template_echo) '_bold_template.nii']);

    % Create cell array of scan names to realign template echo timeseries
    functional_fn = fullfile(preproc_dir, sub, 'func', [previous_steps_prefix sub '_task-' task '_run-' num2str(run) '_echo-' num2str(reference_echo) '_bold.nii']);
    scans = {};
    scans{1} = template_vol;
    for i = 2:N_vol + 1
        scans{i} = [functional_fn ',' num2str(i - 1)];
    end

    % Realign (estimate and reslice) template echo timeseries
    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    realign_estimate_reslice = struct;
    % Data
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.data = {scans'};
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
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = prefix;
    % Run
    spm_jobman('run', realign_estimate_reslice.matlabbatch);
    [d, fn, ext] = fileparts(functional_fn);
    rfunctional_fn = [d filesep prefix fn ext];

    % Access realignment (head motion) parameters
    [d, fn, ext] = fileparts(template_vol);
    HMP_temp_fn = [d filesep 'rp_' fn '.txt'];
    HMP = load(HMP_temp_fn);
    HMP(1, :) = [];
    HMP_fn = fullfile(d, ['rp_' sub '_task-' task '_run-' num2str(run) '.txt']);
    dlmwrite(HMP_fn, HMP, 'delimiter', '\t', 'precision', '%1.7e');
    delete(HMP_temp_fn);

    for e = 1:N_echoes
        if e == reference_echo
            continue
        end
        % Apply rigid body transformation, estimated from HMPs, to each volume in timeseries
        % When applying the transformation, new data is saved to the image header. If applied to the original timeries,
        % this will change data without changing the filename. For this reason, the transformations are applied to a copy of
        % the timeseries, with a prefix 'r'
        functional_fn = fullfile(preproc_dir, sub, 'func', [previous_steps_prefix sub '_task-' task '_run-' num2str(run) '_echo-' num2str(e) '_bold.nii']);
        rfunctional_fn = fullfile(preproc_dir, sub, 'func', [prefix previous_steps_prefix sub '_task-' task '_run-' num2str(run) '_echo-' num2str(e) '_bold.nii']);
        functional_spm = spm_vol(functional_fn);
        functional_img = spm_read_vols(functional_spm);
        rfunctional_spm = functional_spm;
        for i = 1:N_vol
            currentVol = functional_spm(i);
            currentImg = functional_img(:, :, :, i);
            Pm = zeros(12, 1);
            Pm(1:6) = HMP(i, :);
            orig_mat = currentVol.mat;
            rigid_mat = spm_matrix(Pm, 'T*R');
            trans_mat = rigid_mat * orig_mat;
            rfunctional_spm(i).mat = trans_mat;
            rfunctional_spm(i).fname = rfunctional_fn;
            rfunctional_spm(i).private.dat.fname = rfunctional_fn;
            spm_write_vol(rfunctional_spm(i), functional_img(:, :, :, i));
        end

        % Reslice all volumes with reference to template volume
        scans = {};
        scans{1} = template_vol;
        for i = 2:N_vol + 1
            scans{i} = [rfunctional_fn ',' num2str(i - 1)];
        end
        % Reslice echo timeseries
        spm('defaults', 'fmri');
        spm_jobman('initcfg');
        realign_reslice = struct;
        % Data
        realign_reslice.matlabbatch{1}.spm.spatial.realign.write.data = scans';
        % Roptions
        realign_reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.which = [1 0]; % images [2..n]
        realign_reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
        realign_reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
        realign_reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
        realign_reslice.matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = prefix;
        % Run
        spm_jobman('run', realign_reslice.matlabbatch);
        [d, fn, ext] = fileparts(rfunctional_fn);
        rrfunctional_fn = [d filesep prefix fn ext];
        % After reslicing, another prefix 'r' is added. Thus, the original timeseries now has an 'rr' prefix. To
        % keep the naming conventions standard, the filename with 'r' prefix (updated transormation matrices) is deleted
        % and the filename with 'rr' prefix (resliced images) is renamed to have 'r' prefix.
        %% Delete
        delete(rfunctional_fn);
        [status, msg, msgID] = movefile(rrfunctional_fn, rfunctional_fn);
        if status == 0
            disp(msg);
        end
    end
