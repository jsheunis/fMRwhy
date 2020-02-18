function f = displayMaskContour(epi_img, mask, rotate, Nskip)

% Simple function to create a montage / mosaic of multiple slices from a
% single 3D EPI image matrix (typically the mean EPI from a time series)
% with the calculated mask contour as an overlay. This enables quality
% checking of the segmentation and maskin procedure. Needs improvement
% i.t.o RAS/LAS orientation specification and image layout...
%
% INPUT:
% epi_img       - 3D (x,y,z) EPI image matrix
% mask          - 3D (x,y,z) binary mask image matrix
% rotate        - rotate images 90 deg clockwise? yes = 1; no = 0.
% Nskip         - number of slices to skip in between displayed slices
% 
% OUTPUT: 
% output        - structure with montage data
%__________________________________________________________________________
% Copyright (C) Stephan Heunis 2018

[Ni, Nj, Nk] = size(epi_img);
% Rotate image slices if required
if rotate
    for p = 1:Nk
        epi_img(:,:,p) = rot90(epi_img(:,:,p));
    end
end

% Get the boundaries of included mask voxels per slice
bound = cell(Nk,1);
for k = 1:Nk
    bound{k,1} = bwboundaries(squeeze(mask(:, :, k)));
end

% Create figure
f = figure;
im = cell(9,1);
% 1 - For each element in a 3 column 3 row matrix, concatenate multiple epi
% image slices. Plot.
% 2 - For same elements in a 3 column 3 row matrix, concatenate correct
% boundary/contour image slices. Plot as overlay.
for r = 1:3
    for c = 1:3
        el = 3*(r-1)+c;
        subplot(3,3,el);
        im{el,1} = imagesc(epi_img(:, :, Nskip*el));
        colormap gray;
        hold on;
        bmax = numel(bound{Nskip*el, 1});
        for b = 1:bmax
            ax = plot(bound{Nskip*el,1}{b,1}(:,2), bound{Nskip*el,1}{b,1}(:,1), 'r', 'LineWidth', 2.5);
            if rotate
                direction = [0 0 1];
                rotate(ax,direction,90);
            end
        end
        hold off;
    end
end

subplot(3,3,2); title('Brain mask contours overlayed on mean EPI')
