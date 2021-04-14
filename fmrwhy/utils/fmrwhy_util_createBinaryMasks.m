function [GM_img_bin, WM_img_bin, CSF_img_bin, brain_img_bin] = fmrwhy_util_createBinaryMasks(gm_fn, wm_fn, csf_fn, threshold)

    % DESCRIPTION:
    % This function constructs 3D binary masks for GM, WM and CSF based on the
    % relative value of GM/WM/CSF tissue probability maps per voxel. It assumes
    % that the image filename parameters are for images that are in the same
    % space (i.e. they match in size, voxel by voxel). If a threshold is
    % specified (0<= threshold <=1), the images are first thresholded before
    % the binary masks are calculated. For no threshold, specify zero.
    %
    % INPUT:
    % gm_fn         - filename of pre-real-time functional scan
    % wm_fn         - filename of T1-weighted structural scan
    % csf_fn        - filename of T1-weighted structural scan
    % threshold     -
    %
    % OUTPUT:
    % GM_img_bin        - structure with filenames and data
    % WM_img_bin        - structure with filenames and data
    % CSF_img_bin        - structure with filenames and data
    % brain_img_bin        - structure with filenames and data
    % __________________________________________________________________________
    % Copyright (C) 2018 Neu3CA.org
    % Written by Stephan Heunis

    nii = nii_tool('load', gm_fn);
    GM_img = nii.img;
    nii = nii_tool('load', wm_fn);
    WM_img = nii.img;
    nii = nii_tool('load', csf_fn);
    CSF_img = nii.img;

    if threshold ~= 0
        GM_img_thresh = GM_img;
        WM_img_thresh = WM_img;
        CSF_img_thresh = CSF_img;
        GM_img_thresh(GM_img < threshold) = 0;
        WM_img_thresh(WM_img < threshold) = 0;
        CSF_img_thresh(CSF_img < threshold) = 0;
        GM_img_bin = (GM_img_thresh >= WM_img_thresh) & (GM_img_thresh >= CSF_img_thresh) & (GM_img_thresh ~= 0);
        WM_img_bin = (WM_img_thresh > GM_img_thresh) & (WM_img_thresh >= CSF_img_thresh) & (WM_img_thresh ~= 0);
        CSF_img_bin = (CSF_img_thresh > GM_img_thresh) & (CSF_img_thresh > WM_img_thresh) & (CSF_img_thresh ~= 0);
    else
        GM_img_bin = (GM_img >= WM_img) & (GM_img >= CSF_img) & (GM_img ~= 0);
        WM_img_bin = (WM_img > GM_img) & (WM_img >= CSF_img) & (WM_img ~= 0);
        CSF_img_bin = (CSF_img > GM_img) & (CSF_img > WM_img) & (CSF_img ~= 0);
    end

    brain_img_bin = GM_img_bin | WM_img_bin | CSF_img_bin;
