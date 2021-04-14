function fmrwhy_util_saveNiftiFrom4D(fn4D, fn3D, idx)

    % Create a single 3D nifti file from a 4D nifti file, given a specific image index.fn3D
    % OR alternatively, split 4D nifti file into a full set of 3D nifti images

    % Load template 4D nii
    template_nii = nii_tool('load', fn4D);

    if nargin < 3
        % Convert 4d nifti to set of 3d niftis
        nii_tool('save', template_nii, fn3D, true);
    else
        % Create new nifti with hdr of template nii
        new_nii = struct;
        new_nii.hdr = template_nii.hdr;
        new_nii.img = template_nii.img(:, :, :, idx);
        new_nii.hdr.dim(1) = 3; % set number of dimensions of image ==> 3D (not 4D)
        new_nii.hdr.dim(5) = 1; % set the value of 4th dimension (time) to 1
        new_nii.hdr.aux_file = '';
        new_nii.hdr.file_name = fn3D;
        nii_tool('save', new_nii, fn3D);
    end
