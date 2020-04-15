function fmrwhy_util_saveNifti(new_fn, new_img, template_fn)

% Load template nii
template_nii = nii_tool('load', template_fn);

% Create new nifti with hdr of template nii
new_nii = struct;
new_nii.hdr = template_nii.hdr;
new_nii.img = new_img;
new_nii.hdr.aux_file = '';
new_nii.hdr.file_name = new_fn;
nii_tool('save', new_nii, new_fn);