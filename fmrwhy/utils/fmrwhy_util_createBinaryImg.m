function img_bin = fmrwhy_util_createBinaryImg(img, threshold)

    % DESCRIPTION:
    % This function constructs 3D binary mask for the input image. If a threshold is
    % specified (0<= threshold <=1), the images are first thresholded before
    % the binary masks are calculated. For no threshold, specify zero.
    %
    % __________________________________________________________________________
    % Copyright (C) 2018

    if threshold ~= 0
        img_thresh = img;
        img_thresh(img < threshold) = 0;
        img_bin = img_thresh ~= 0;
    else
        img_bin = img ~= 0;
    end
