function fmrwhy_util_applyTransform(toTransform_fn, motion_params, template_fn, saveAs_fn)

[Ni, Nj] = size(motion_params);
functional_spm = spm_vol(toTransform_fn);
N_vol = numel(functional_spm);

if Ni ~= N_vol
    disp(['Error: dimension mismatch between head movement parameters (' num2str(Ni) 'points) and functional timeseries (' num2str(N_vol) 'points)' ]);
    return;
end

functional_img = spm_read_vols(functional_spm);
tfunctional_spm = functional_spm;

% For each volume, get transformation matrix and apply to .mat
for i = 1:N_vol
    currentVol = functional_spm(i);
    currentImg = functional_img(:,:,:,i);
    Pm = zeros(12,1);
    Pm(1:6) = motion_params(i, :);
    orig_mat = currentVol.mat;
    rigid_mat = spm_matrix(Pm, 'T*R');
    trans_mat = rigid_mat * orig_mat;
    tfunctional_spm(i).mat = trans_mat;
    tfunctional_spm(i).fname = saveAs_fn;
    tfunctional_spm(i).private.dat.fname = saveAs_fn;
    spm_write_vol(tfunctional_spm(i),functional_img(:,:,:,i));
end

% Reslice all volumes with reference to template volume
reslice_scans = {};
reslice_scans{1} = [template_fn ',1'];
for i = 2:N_vol+1
    reslice_scans{i} = [saveAs_fn ',' num2str(i-1)];
end
fmrwhy_batch_realignResl(reslice_scans)

[d, fn, ext] = fileparts(saveAs_fn);
rsaveAs_fn = [d filesep 'rr' fn ext];
% After reslicing, a prefix 'rr' is added.
% First, the saveAS_fn filename (with updated transormation matrices) is deleted
% Then, the filename with 'rr' prefix (resliced images) is renamed to saveAS_fn
delete(saveAs_fn)
[status, msg, msgID] = movefile(rsaveAs_fn, saveAs_fn);
if status == 0
    disp(msg)
end



