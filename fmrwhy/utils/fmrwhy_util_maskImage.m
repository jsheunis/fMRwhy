function out_img = fmrwhy_util_maskImage(img, mask)

    % img = 3D array
    % mask = binary 3D array

    in_img = img(:);
    out_img = zeros(size(in_img));
    I_mask = find(mask(:));
    out_img(I_mask) = in_img(I_mask);
    out_img = reshape(out_img, size(img));
