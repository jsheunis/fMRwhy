function f = getSPM1stLevel(stats_dir)

f = struct;
SPM_fn = [stats_dir filesep 'SPM.mat'];
SPM_mat = load(SPM_fn);

f.X_design = SPM_mat.SPM.xX.X;
f.Nx = SPM_mat.SPM.xVol.DIM(1);
f.Ny = SPM_mat.SPM.xVol.DIM(2);
f.Nz = SPM_mat.SPM.xVol.DIM(3);

[f.Nt, f.Nregr] = size(SPM_mat.SPM.xX.X);
f.beta_img = {};
f.beta = zeros(f.Nregr, f.Nx*f.Ny*f.Nz);

for i = 1:f.Nregr
    beta_fn = [stats_dir filesep 'beta_' sprintf('%04d',i) '.nii'];
    f.beta_img{i,1} = spm_read_vols(spm_vol(beta_fn));
    f.beta(i, :) = reshape(f.beta_img{i,1}, 1, f.Nx*f.Ny*f.Nz);
end

mask_fn = [stats_dir filesep 'mask.nii'];
f.mask_img = spm_read_vols(spm_vol(mask_fn));
f.mask = reshape(f.mask_img, f.Nx*f.Ny*f.Nz, 1);
f.I_mask = find(f.mask);


