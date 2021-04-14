function fmrwhy_util_applyTransform(toTransform_fn, motion_params, template_fn, saveAs_fn)

    [Ni, Nj] = size(motion_params);
    functional_spm = spm_vol(toTransform_fn);
    N_vol = numel(functional_spm);

    if Ni ~= N_vol
        disp(['Error: dimension mismatch between head movement parameters (' num2str(Ni) 'points) and functional timeseries (' num2str(N_vol) 'points)']);
        return
    end

    functional_img = spm_read_vols(functional_spm);
    tfunctional_spm = functional_spm;

    flagsSpmReslice = struct('quality', .9, 'fwhm', 5, 'sep', 4, 'interp', 4, 'wrap', [0 0 0], 'mask', 1, 'mean', 0, 'which', 2);

    funcref_spm = spm_vol(template_fn);
    funcref_3D = spm_read_vols(funcref_spm);
    R = struct;
    R(1, 1).mat = funcref_spm.mat;
    R(1, 1).dim = funcref_spm.dim;
    R(1, 1).Vol = funcref_3D;
    reslVol = cell(N_vol, 1);

    % For each volume, get transformation matrix and apply to .mat
    for i = 1:N_vol
        currentVol = functional_spm(i);
        currentImg = functional_img(:, :, :, i);
        Pm = zeros(12, 1);
        Pm(1:6) = motion_params(i, :);
        orig_mat = currentVol.mat;
        rigid_mat = spm_matrix(Pm, 'T*R');
        trans_mat = rigid_mat * orig_mat;

        R(2, 1).mat = trans_mat;
        R(2, 1).dim = currentVol.dim;
        R(2, 1).Vol = currentImg;

        reslVol{i} = fmrwhy_realtime_reslice(R, flagsSpmReslice, currentVol);

        tfunctional_spm(i).fname = saveAs_fn;
        tfunctional_spm(i).private.dat.fname = saveAs_fn;
        spm_write_vol(tfunctional_spm(i), reslVol{i});

        % TODO: replace this with dicm2nii equivalent (ask on github repo)
    end
