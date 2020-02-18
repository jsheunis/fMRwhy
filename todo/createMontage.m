function montage = createMontage(img, columns, rotate, str, clrmp)

% Simple function to create a montage / mosaic of multiple slices from a single 3D
% image matrix. Needs improvement i.t.o RAS/LAS orientation specification
% and image layout...
%
% INPUT:
% img           - 3D (x,y,z) image matrix
% columns       - number of columns in montage (rows are calculated accordingly)
% rotate        - rotate images 90 deg clockwise? yes = 1; no = 0.
% str           - figure title
% clrmp         - figure colormap

% 
% OUTPUT: 
% output        - structure with montage data
%__________________________________________________________________________
% Copyright (C) Stephan Heunis 2018


montage = struct;
[Ni, Nj, Nk] = size(img);
% Rotate image slices if required
if rotate
    for p = 1:Nk
        img(:,:,p) = rot90(img(:,:,p));
    end
end

% Determine amount of rows and filler slices
rows = floor(Nk/columns);
fill = mod(Nk, columns);
filler = zeros(Ni, Nj);
if fill == 0
    N_fill = 0;
else
    N_fill = columns - mod(Nk, columns);
end
montage.rows = rows;
montage.columns = columns;
montage.N_fill = N_fill;

parts = {};
% 1 - Concatenate slices together horizontally, per row (except last).
% 2 - Concatenate rows together vertically
for i = 1:rows
    for j = 1:columns
        if j ==1
            parts{i} = img(:,:,(columns*(i-1)+j));
        else
            parts{i} = cat(2, parts{i}, img(:,:,(columns*(i-1)+j)));
        end
    end
    if i ==1
        whole = parts{i};
    else
        whole = cat(1, whole, parts{i});
    end
end

% 1 - Concatenate filler slices to last row, if required.
% 2 - Concatenate last row to whole matrix, if required.
if N_fill ~= 0
    % last row
    last_parts = img(:,:,(rows*columns+1));
    for k = (rows*columns+2):Nk
        last_parts = cat(2, last_parts, img(:,:,k));
    end
    for m = 1:N_fill
        last_parts = cat(2, last_parts, filler);
    end
    montage.image = cat(1, whole, last_parts);
else
    montage.image = whole;
end

% Create figure
f = figure;imagesc(montage.image); colormap(clrmp); colorbar;
title(str);
montage.f = f;





