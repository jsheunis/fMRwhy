function output = calculateQCmeasures(functional4D_fn, structural_fn, fwhm, spm_dir, out_dir, subject)

% NOTE: THIS TOOL IS UNDER CONTINUOUS DEVELOPMENT AND HAS NOT BEEN SUFFCIENTLY TESTED

% Function to calculate multiple quality control measures for an fMRI time
% series using SPM12 and Matlab. The goal is to create a tool, for use by
% fMRI technicians/researchers/clinicians familiar with SPM12 and Matlab,
% that allows the calculation of measures that can assist in diagnosing
% quality issues in fMRI data. It is not intended for use as ground truth
% for quality diagnosis, as these measures are known to be relative and to
% vary based on scanner site, acquisition time, data format, and more.
% However, it can provide insight into possible data quality issues
% originating from the scanner or single subject.

% Function steps include preprocessing the acquired fMRI and structural data
% (coregistering structural image to first functional image, segmenting the
% coregistered structural image into tissue types, and reslicing the
% segments to the functional resolution image grid), and calling functions
% for a standard set of QC measures based on various sources (mainly
% including the Quality Assessment Protocol (QAP) from the Preprocessed
% Connectomes Project (PCP) and MRIQC. This function makes use of SPM12
% functions and batch routines. If SPM12 batch parameters are not
% explicitly set, defaults are assumed. 
% 
% CURRENT QC MEASURES (27/06/2018):
% - Temporal signal to noise ratio (tSNR) (3D image)
% - Mean tSNR of brain, grey matter, white matter, and CSF
% - Z-score timeseries and mean
% - Standard deviation (3D image)
% - Framewise displacement (FD) timeseries, total and mean (mm)
% - GCOR = global correlation ....
% - Non-standardized differential variance (DVARS) timeseries (a.u.)
% - The Plot (GM, WM and CSF voxel intensities over time)
%
% INPUT:
% funcional4D_fn     - filename of 4D functional timeseries (.nii only)
% structural_fn      - filename of T1-weighted structural scan (.nii only)
% fwhm               - kernel size for smoothing operations (mm)
% spm_dir            - SPM12 directory
% out_dir            - output directory (for figures and html logfile)
% subject            - subject name/code
%
% OUTPUT:
% Matlab figures:   - timeseries plots (FD, DVARS, Zscore, The Plot)
%                   - montage images (tSNR, tSNR brain, stddev, mean EPI)
%                   - coregistration contour plot
% HTML logfile:     - HTML logfile with measures and figures and metadata
%
% SOURCES:
% PCP-QAP:          - http://preprocessed-connectomes-project.org/quality-assessment-protocol/index.html
% MRIQC:            - https://mriqc.readthedocs.io/en/latest/
% The Plot:         - https://www.sciencedirect.com/science/article/pii/S1053811916303871?via%3Dihub
%__________________________________________________________________________
% Copyright (C) Stephan Heunis 2018


% User defined variables
% -------------------------------------------------------------------------
intensity_scale = [-6 6]; % scaling for plot image intensity, see what works
FD_threshold = 0.5; % mm
% -------------------------------------------------------------------------
if ~exist(out_dir, 'dir')
    mkdir(out_dir)
end
cd(out_dir)

output = struct;

% Get image information
func_spm = spm_vol(functional4D_fn);
tsize = size(func_spm);
Nt = tsize(1);
Ni= func_spm(1).dim(1);
Nj= func_spm(1).dim(2);
Nk= func_spm(1).dim(3);

% Preprocess structural and functional data.
[d, f, e] = fileparts(structural_fn);
if exist([d filesep 'rc1' f e], 'file')
    % If a resliced grey matter image is present in the specified data
    % directory, assume that all preprocessing has been completed by
    % standardPreproc, and declare appropriate variables.
    disp('Preproc done!')
    preproc_data = struct;
    preproc_data.forward_transformation = [d filesep 'y_' f e];
    preproc_data.inverse_transformation = [d filesep 'iy_' f e];
    preproc_data.gm_fn = [d filesep 'c1' f e];
    preproc_data.wm_fn = [d filesep 'c2' f e];
    preproc_data.csf_fn = [d filesep 'c3' f e];
    preproc_data.bone_fn = [d filesep 'c4' f e];
    preproc_data.soft_fn = [d filesep 'c5' f e];
    preproc_data.air_fn = [d filesep 'c6' f e];
    preproc_data.rstructural_fn = [d filesep 'r' f e];
    preproc_data.rgm_fn = [d filesep 'rc1' f e];
    preproc_data.rwm_fn = [d filesep 'rc2' f e];
    preproc_data.rcsf_fn = [d filesep 'rc3' f e];
    preproc_data.rbone_fn = [d filesep 'rc4' f e];
    preproc_data.rsoft_fn = [d filesep 'rc5' f e];
    preproc_data.rair_fn = [d filesep 'rc6' f e];
    
    [d, f, e] = fileparts(functional4D_fn);
    preproc_data.rfunctional_fn = [d filesep 'r' f e];
    preproc_data.srfunctional_fn = [d filesep 'sr' f e];
    preproc_data.sfunctional_fn = [d filesep 's' f e];
    preproc_data.mp_fn = [d filesep 'rp_' f '.txt'];
    preproc_data.MP = load(preproc_data.mp_fn);
else
    % If a resliced grey matter image is NOT present in the specified data
    % directory, run standardPreproc
    preproc_data = standardPreproc(functional4D_fn, structural_fn, fwhm, spm_dir);
end

% Calculate brain mask matrices for GM, WM, CSF, and all together
mask_threshold = 0.5;
[GM_img_bin, WM_img_bin, CSF_img_bin] = createBinarySegments(preproc_data.rgm_fn, preproc_data.rwm_fn, preproc_data.rcsf_fn, mask_threshold);
I_GM = find(GM_img_bin);
I_WM = find(WM_img_bin);
I_CSF = find(CSF_img_bin);
mask_reshaped = GM_img_bin | WM_img_bin | CSF_img_bin;
I_mask = find(mask_reshaped);
Nmaskvox = numel(I_mask);

% Detrend 4D time series
output = detrend4D(functional4D_fn);
F4D_detrended = output.F4D_detrended;
F2D_detrended = output.F_2D_detrended;

% Statistical measures
F2D_mean = mean(F2D_detrended, 2);
F2D_stddev = std(F2D_detrended, 0, 2);
F2D_zstat = (F2D_detrended - F2D_mean)./F2D_stddev;
F2D_zstat(isnan(F2D_zstat))=0;
F2D_zstat_mean = mean(abs(F2D_zstat),1);
Zstat_mean = mean(F2D_zstat_mean);
F2D_var = var(F2D_detrended,0,2);
F2D_psc = 100*(F2D_detrended./repmat(F2D_mean, 1, Nt)) - 100;
F2D_psc(isnan(F2D_psc))=0;

% Framewise displacement
r = 80; % mm
FD_measures = calculateFD(preproc_data.MP, r, FD_threshold);

% DVARS
F2D_diff = [zeros(1, Ni*Nj*Nk); diff(F2D_detrended')]';
DVARS = var(F2D_diff);

% GCOR
% Steps according to https://doi.org/10.1089/brain.2013.0156:
% (1)?De-mean each voxel's time series and scale it by its Eucledian norm
% (2)?Average scaled time series over the whole brain mask
% (3)?GCOR is the length (L2 norm) of this averaged series
F2D_demeaned = F2D_detrended - F2D_mean;
F2D_norms = sqrt(sum(F2D_demeaned.^2,2));
F2D_scaled = F2D_demeaned./F2D_norms;
F2D_ave_timeseries = mean(F2D_scaled(I_mask,:), 1);
GCOR = sum(F2D_ave_timeseries.^2,2);


% The Plot (with FD, DVARS, and mean Zscore per volume)
GM_img = F2D_psc(I_GM, :);
WM_img = F2D_psc(I_WM, :);
CSF_img = F2D_psc(I_CSF, :);
all_img = [GM_img; WM_img; CSF_img];
line1_pos = numel(I_GM);
line2_pos = numel(I_GM) + numel(I_WM);
tf = figure;
fontsizeL = 14;
fontsizeM = 11;
ax1 = subplot(7,1,4:7);
imagesc(ax1, all_img); colormap(gray); caxis(intensity_scale);
title(ax1, 'thePlotSpm','fontsize',fontsizeL)
ylabel(ax1, 'Voxels','fontsize',fontsizeM)
xlabel(ax1, 'fMRI volumes','fontsize',fontsizeM)
hold on; line([1 Nt],[line1_pos line1_pos],  'Color', 'b', 'LineWidth', 2 )
line([1 Nt],[line2_pos line2_pos],  'Color', 'r', 'LineWidth', 2 )
hold off;
ax2 = subplot(7,1,1);
plot(ax2, FD_measures.FD, 'LineWidth', 2); grid;
set(ax2,'Xticklabel',[]);
title(ax2, 'FD','fontsize',fontsizeL)
ylabel(ax2, 'mm','fontsize',fontsizeM)
ax3 = subplot(7,1,2);
plot(ax3, DVARS, 'LineWidth', 2); grid;
set(ax3,'Xticklabel',[]);
title(ax3, 'DVARS','fontsize',fontsizeL)
ylabel(ax3, 'a.u.','fontsize',fontsizeM)
ax4 = subplot(7,1,3);
plot(ax4, F2D_zstat_mean, 'LineWidth', 2); grid;
set(ax4,'Xticklabel',[]);
title(ax4, 'Z-score','fontsize',fontsizeL)
ylabel(ax4, 'a.u.','fontsize',fontsizeM)
print(tf, 'timeseries_summary', '-dpng')

output.all_img = all_img;
output.F2D_psc = F2D_psc;
output.I_GM = I_GM;
output.I_WM = I_WM;
output.I_CSF = I_CSF;
output.I_mask = I_mask;
output.Nmaskvox = Nmaskvox;
output.Ni = Ni;
output.Nj = Nj;
output.Nk = Nk;
output.Nt = Nt;

% tSNR
tSNR_2D = F2D_mean./F2D_stddev;
tSNR_brain = mean(tSNR_2D(I_mask));
tSNR_GM = mean(tSNR_2D(I_GM));
tSNR_WM = mean(tSNR_2D(I_WM));
tSNR_CSF = mean(tSNR_2D(I_CSF));

% Display metrics in command window
disp(['Number of volumes classified as outliers based on FD>=' num2str(FD_threshold) 'mm: ' num2str(numel(FD_measures.FD_outliers_ind))])
disp(['Total FD: ' num2str(FD_measures.FD_sum)])
disp(['Mean FD: ' num2str(FD_measures.FD_mean)])
disp(['Mean Zscore: ' num2str(Zstat_mean)])
disp(['GCOR: ' num2str(GCOR)])
disp(['tSNR (brain): ' num2str(tSNR_brain)])
disp(['tSNR (GM): ' num2str(tSNR_GM)])
disp(['tSNR (WM): ' num2str(tSNR_WM)])
disp(['tSNR (CSF): ' num2str(tSNR_CSF)])

% Prepare 3D and 4D images
mask_3D = reshape(mask_reshaped, Ni, Nj, Nk);
tSNR_3D = reshape(tSNR_2D, Ni, Nj, Nk);
F3D_mean = reshape(F2D_mean, Ni, Nj, Nk);
F3D_var = reshape(F2D_var, Ni, Nj, Nk);
F3D_stddev = reshape(F2D_stddev, Ni, Nj, Nk);
tSNR_2D_masked = zeros(Ni*Nj*Nk, 1);
tSNR_2D_masked(I_mask, :) = tSNR_2D(I_mask, :);
tSNR_3D_masked = reshape(tSNR_2D_masked, Ni, Nj, Nk);

% Create montages of 3D images
montage2 = createMontage(F3D_mean, 5, 1, 'Mean EPI (whole image)', 'gray');
print(montage2.f, 'mean_epi', '-dpng')
montage3 = createMontage(F3D_stddev, 5, 1, 'Standard deviation (whole image)', 'parula');
print(montage3.f, 'stddev_epi', '-dpng')
montage1 = createMontage(tSNR_3D, 5, 1, 'tSNR (whole image)', 'hot');
print(montage1.f, 'tsnr_whole', '-dpng')
montage4 = createMontage(tSNR_3D_masked, 5, 1, 'tSNR (brain)', 'hot');
print(montage4.f, 'tsnr_brain', '-dpng')
figmask = displayMaskContour(F3D_mean, mask_3D, 0, 3);
print(figmask, 'mask_contour', '-dpng')


% Save tSNR_3D image - method 1 (spmup Pernet tsnr method)
V = func_spm(1);
V.fname = [out_dir filesep 'tSNR_' subject '_meth1.nii'];
V.private.dat.fname  = [out_dir filesep 'tSNR_' subject '_meth1.nii'];
V.descrip = 'tSNR image - see spmup_temporalSNR';
V.private.dat.dim = V(1).private.dat.dim(1:3);
V.n = [1 1];
spm_write_vol(V,tSNR_3D);

% Save tSNR_3D image - method 2 (mixture of spmup Pernet tsnr method, and info from here: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;12fa60a.1205)
Va = func_spm(1);
Va.fname = [out_dir filesep 'tSNR_' subject '_meth2.nii'];
Va.private.dat.fname  = [out_dir filesep 'tSNR_' subject '_meth2.nii'];
Va.descrip = 'tSNR image - see spmup_temporalSNR';
Va.private.dat.dim = Va(1).private.dat.dim(1:3);
Va.n = [1 1];
Va.pinfo(1) = 1; % https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;12fa60a.1205
spm_write_vol(Va,tSNR_3D);

% Save tSNR_3D image - method 3
Vb = func_spm(1);
new_nii = spm_create_vol(Vb);
new_nii.fname = [out_dir filesep 'tSNR_' subject '_meth3.nii'];
new_img = spm_write_vol(new_nii, tSNR_3D);



% Create HTML log file
dt = datetime('now');
[Y,MO,D,H,MI,S] = datevec(dt);
dt_str = [num2str(Y) num2str(MO) num2str(D) num2str(H) num2str(MI) num2str(round(S))];
t = datestr(dt);

log_name = [subject '_' dt_str '.html'];
fid = fopen(log_name,'a');
fprintf(fid, '<H2>Log</H2>');
fprintf(fid, ['\n<BR>Subject:  ' subject]);

fprintf(fid, ['\n<BR>Date/time:  ' t]);

fprintf(fid, '<H2>Imaging info</H2>');
fprintf(fid, ['\nVolumes:  ' num2str(Nt)]);
fprintf(fid, ['\n<BR>Voxels (x,y,z):  ' num2str(Ni) ', ' num2str(Nj) ', ' num2str(Nk)]);

fprintf(fid, '<H2>Timeseries summary</H2>');
fprintf(fid, '\n<TABLE><TR><TD><img src="timeseries_summary.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );

fprintf(fid, '<H2>QC Metrics</H2>');
fprintf(fid, ['\nFD threshold (mm):  ' num2str(FD_threshold)]);
fprintf(fid, ['\n<BR>FD outliers:  ' num2str(numel(FD_measures.FD_outliers_ind))]);
fprintf(fid, ['\n<BR>Total FD:  ' num2str(FD_measures.FD_sum)]);
fprintf(fid, ['\n<BR>Mean FD:  ' num2str(FD_measures.FD_mean)]);
fprintf(fid, ['\n<BR>Mean Zscore:  ' num2str(Zstat_mean)]);
fprintf(fid, ['\n<BR>GCOR:  ' num2str(GCOR)]);
fprintf(fid, ['\n<BR>tSNR (brain):  ' num2str(tSNR_brain)]);
fprintf(fid, ['\n<BR>tSNR (GM):  ' num2str(tSNR_GM)]);
fprintf(fid, ['\n<BR>tSNR (WM):  ' num2str(tSNR_WM)]);
fprintf(fid, ['\n<BR>tSNR (CSF):  ' num2str(tSNR_CSF)]);

fprintf(fid, '<H2>QC brain images</H2>');
fprintf(fid, '\n<TABLE><TR><TD><img src="mean_epi.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );
fprintf(fid, '\n<TABLE><TR><TD><img src="stddev_epi.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );
fprintf(fid, '\n<TABLE><TR><TD><img src="tsnr_whole.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );
fprintf(fid, '\n<TABLE><TR><TD><img src="tsnr_brain.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );
fprintf(fid, '\n<TABLE><TR><TD><img src="mask_contour.png" alt="no picture" width=700 height=700></TD></TR></TABLE>' );
fclose(fid);


























