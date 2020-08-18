function output = fmrwhy_util_getBinaryOverlap(binary_fns, saveAs_overlap_fn, background_fn, saveAs_montage_fn)

% DESCRIPTION:
%
% INPUT:
%
% OUTPUT:
%__________________________________________________________________________
% Copyright (C)
% Written by Stephan Heunis

N_imgs = numel(binary_fns);
%N_imgs = numel(binary_imgs);
if N_imgs < 2
    disp('Error: overlap can not be determined for fewer than 2 images')
    return;
end

output = struct;
output.size_imgs = {};
output.binary_imgs = {};
for i = 1:N_imgs
    [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(binary_fns{i});
    output.binary_imgs{i} = fmrwhy_util_createBinaryImg(p.nii.img, 0.1); % Todo: decide if this binarization should be done outside of this function, in which case images (not fns) need to be passed as parameters
    output.size_imgs{i} = sum(output.binary_imgs{i}(:));
    if i == 1
        output.overlap = output.binary_imgs{i};
    else
        output.overlap = output.overlap & output.binary_imgs{i};
    end
end

% Determine overlap size
output.size_overlap = sum(output.overlap(:));

% Dice coefficient if 2 images
if N_imgs == 2
    output.dice_coeff = 2*output.size_overlap/(output.size_imgs{1} + output.size_imgs{2});
end

% Fractions: overlap vs images
output.fraction = {}
for i = 1:N_imgs
    output.fraction{i} = output.size_overlap / output.size_imgs{i};
end

%
%no_scaling = 1;
%str = consess{k}.tcon.name;
%
%if nargin == 2
%
%saveAs_overlap_fn = fullfile(run_dir_stats, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-' str 'Overlaps' options.roi.(task).desc{j} '_roi.nii']);
%fmrwhy_util_saveNifti(saveAs_overlap_fn, double(overlap_img), tmap_clusters_fn, no_scaling);
%[poverlap, ~, ~, ~] = fmrwhy_util_readOrientNifti(saveAs_overlap_fn);
%
%plot_overlap_img = double(poverlap.nii.img);
%saveAs_fn = strrep(saveAs_overlap_fn, '.nii', '.png');
%roi_rgbcolors = [148, 239, 255];
%
%overlapmontage = fmrwhy_util_createStatsOverlayMontage(background_img, [], plot_overlap_img, 9, 1, '', 'gray', 'off', 'max', [], [], roi_rgbcolors, false, saveAs_fn);