function fmrwhy_util_saveNifti(new_fn, new_img, template_fn, no_scaling)

% no_scaling: if 1 or true, do not scale image values according to template values
% useful when saving quantitative images

if nargin < 4
    no_scaling = 0;
end

% Load template nii
template_nii = nii_tool('load', template_fn);

% Create new nifti with hdr of template nii
new_nii = struct;
new_nii.hdr = template_nii.hdr;
new_nii.img = new_img;
new_nii.hdr.aux_file = '';
new_nii.hdr.file_name = new_fn;
if no_scaling
    new_nii.hdr.scl_slope = [];
    new_nii.hdr.scl_inter = [];
end


nii_tool('save', new_nii, new_fn);