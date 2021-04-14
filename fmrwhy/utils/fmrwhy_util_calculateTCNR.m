function output = fmrwhy_util_calculateTCNR(functional_fn, contrast_vector, mask_fn, saveAs_fn, template_fn)
    % Function to calculate tCNR from (masked) fMRI timeseries data
    % contrast_vector:  binary vector of same length as number of functional time series volumes.
    %                   ones for ON condition, zeros for OFF condition

    % INPUT:

    % INPUT:

    % OUTPUT:
    output = struct;

    % Load timeseries data
    nii = nii_tool('load', functional_fn);
    data_4D = double(nii.img); % [Ni x Nj x Nk x Nt]
    [Ni, Nj, Nk, Nt] = size(data_4D);
    data_2D = reshape(data_4D, Ni * Nj * Nk, Nt); % [voxels, time]
    data_2D = data_2D'; % [time, voxels]

    if Nt ~= numel(contrast_vector)
        disp('Error: contrast_vector must be same length as the functional time series');
        return
    end

    ON_volumes = contrast_vector;
    OFF_volumes = ~contrast_vector;

    % Then get the mask if fn supplied
    if mask_fn ~= 0
        nii_mask = nii_tool('load', mask_fn);
        mask_3D = double(nii_mask.img); % [Ni x Nj x Nk]
        mask_2D = reshape(mask_3D, Ni * Nj * Nk, 1); % [voxels, 1]
        I_mask = find(mask_2D); % [voxels, 1]
        I_mask = I_mask'; % [1, voxels]
    end

    % Remove linear and quadratic trends from data, per voxel
    data_2D_detrended = fmrwhy_util_detrend(data_2D, 2); % [time, voxels]
    output.data_2D_detrended = data_2D_detrended;

    % Calculate means and difference (ignore nan values)
    data_2D_mean_ON = nanmean(data_2D_detrended(ON_volumes)); % [1, voxels]
    data_2D_mean_OFF = nanmean(data_2D_detrended(ON_volumes)); % [1, voxels]
    data_2D_mean = data_2D_mean_ON - data_2D_mean_OFF;
    data_3D_mean = reshape(data_2D_mean, Ni, Nj, Nk);
    output.data_3D_mean = data_3D_mean;

    % Calculate standard deviation
    data_2D_stddev = std(data_2D_detrended); % [1, voxels]
    data_3D_stddev = reshape(data_2D_stddev, Ni, Nj, Nk);
    output.data_3D_stddev = data_3D_stddev;

    % Calculate tCNR
    data_2D_tcnr = data_2D_mean ./ data_2D_stddev; % [1, voxels]
    % data_2D_tcnr(data_2D_tcnr<0) = 0; % TODO: double check if this thresholding is fine to do, or give adequate reasoning
    data_3D_tcnr = reshape(data_2D_tcnr, Ni, Nj, Nk);
    output.data_3D_tcnr = data_3D_tcnr;

    % Mask if required
    if mask_fn ~= 0
        data_2D_tcnr_masked = zeros(size(data_2D_tcnr)); % [1, voxels]
        data_2D_tcnr_masked(I_mask) = data_2D_tcnr(I_mask); % [1, voxels]
        data_3D_tcnr_masked = reshape(data_2D_tcnr_masked, Ni, Nj, Nk);
        output.data_3D_tcnr_masked = data_3D_tcnr_masked;
    end

    % Save to file, if required
    if saveAs_fn ~= 0
        no_scaling = 1;
        if mask_fn ~= 0
            fmrwhy_util_saveNifti(saveAs_fn, data_3D_tcnr_masked, template_fn, no_scaling);
        else
            fmrwhy_util_saveNifti(saveAs_fn, data_3D_tcnr, template_fn, no_scaling);
        end
    end
