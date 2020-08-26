function output = fmrwhy_util_getBinaryOverlap(binary_fns, saveAs_overlap_fn, template_fn, saveAs_montage_fn, background_fn)

% DESCRIPTION:
%
% INPUT:
%
% OUTPUT:
%__________________________________________________________________________
% Copyright (C)
% Written by Stephan Heunis

% Get amount of filenames to use in overlap calc
N_imgs = numel(binary_fns);
if N_imgs < 2
    disp('Error: overlap can not be determined for fewer than 2 images')
    return;
end

output = struct;
output.size_imgs = {};
output.binary_imgs = {};

% For each image:
for i = 1:N_imgs
    % Get image data
    nii = nii_tool('load', binary_fns{i});
    % Make sure the image is in binary format
    output.binary_imgs{i} = fmrwhy_util_createBinaryImg(double(nii.img), 0.1);
    % Count the voxels in the image that are part of the ROI, i.e. voxel value = 1 or TRUE
    output.size_imgs{i} = sum(output.binary_imgs{i}(:));
    % Get logical AND of new image with previous images
    if i == 1
        output.overlap_img = output.binary_imgs{i};
    else
        output.overlap_img = output.overlap_img & output.binary_imgs{i};
    end
end

% Determine overlap voxel count
output.size_overlap = sum(output.overlap_img(:));

% Calculate dice coefficient if 2 images were passed
if N_imgs == 2
    output.dice_coeff = 2*output.size_overlap/(output.size_imgs{1} + output.size_imgs{2});
end

% For each image: calculate the fraction of the overlap vs the image roi size
output.fraction = {}
for i = 1:N_imgs
    output.fraction{i} = output.size_overlap / output.size_imgs{i};
end




% Save overlap image as nifti
if nargin > 1
    no_scaling = 1;
    fmrwhy_util_saveNifti(saveAs_overlap_fn, double(output.overlap_img), template_fn, no_scaling);
end

% Create overlap image montage
if nargin > 3
    [poverlap, ~, ~, ~] = fmrwhy_util_readOrientNifti(saveAs_overlap_fn);
    plot_overlap_img = double(poverlap.nii.img);
    [pbackground, ~, ~, ~] = fmrwhy_util_readOrientNifti(background_fn);
    background_img = double(pbackground.nii.img);
    roi_rgbcolors = [148, 239, 255];
    overlapmontage = fmrwhy_util_createStatsOverlayMontage(background_img, [], plot_overlap_img, 9, 1, '', 'gray', 'off', 'max', [], [], roi_rgbcolors, false, saveAs_montage_fn);
end