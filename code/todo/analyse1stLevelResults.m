%% RTME-FMRI-PROCESSING

%% DEFINITIONS AND VARIABLE DECLARATIONS

data_dir = '/Volumes/Stephan_WD/Kempenhaeghe_data/29052017/FT_Iris_real';
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
addpath(spm_dir);
subj = 'sub-01';
s_fn = [data_dir filesep subj filesep 'anat' filesep subj '_T1w.nii'];
f_fn = [data_dir filesep subj filesep 'func' filesep subj '_task-fingertapping.nii'];
stats_dir = [data_dir filesep subj filesep 'stats'];
fwhm = 4;  % mm

%% PREPROCESSING
[d, f, e] = fileparts(s_fn);
[d1, f1, e1] = fileparts(f_fn);
preproc_data = struct;
% Structural filenames
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
% Functional filenames
preproc_data.rfunctional_fn = [d1 filesep 'r' f1 e1];
preproc_data.srfunctional_fn = [d1 filesep 'sr' f1 e1];
preproc_data.mp_fn = [d1 filesep 'rp_' f1 '.txt'];
preproc_data.MP = load(preproc_data.mp_fn);

%% 1ST LEVEL STATS

load([stats_dir filesep 'SPM.mat']);
[Ntt, Nregr] = size(SPM.xX.X);
sess_params = struct;
sess_params.timing_units = 'scans';
sess_params.timing_RT = 3;
sess_params.cond_name = 'Fingertapping';
sess_params.cond_onset = [11; 31; 51; 71; 91; 111; 131; 151];
sess_params.cond_duration = [10; 10; 10; 10; 10; 10; 10; 10];
contrast_params = struct;
contrast_params.weights = zeros(1, Nregr); 
contrast_params.weights(1) = 1;
contrast_params.name = 'Picture viewing';

spm_runResults_jsh(stats_dir)
pause(5);
save([stats_dir filesep 'xSPM.mat'], 'xSPM')
pause(1);

convolved_task_design = SPM.xX.X(:,1);
drift_regressors = SPM.xX.K.X0;
%%
X_design = [convolved_task_design drift_regressors ones(Ntt,1)];

%% Masking
[GM_img_bin, WM_img_bin, CSF_img_bin] = create_binary_segments(preproc_data.rgm_fn, preproc_data.rwm_fn, preproc_data.rcsf_fn, 0.1);
I_GM = find(GM_img_bin);
I_WM = find(WM_img_bin);
I_CSF = find(CSF_img_bin);
I_total_test = numel(I_GM)+numel(I_WM)+numel(I_CSF);
mask_reshaped = GM_img_bin | WM_img_bin | CSF_img_bin;
I_mask = find(mask_reshaped);
Nmaskvox = numel(I_mask);
Nvox = numel(GM_img_bin);
[Ni, Nj, Nk] = size(GM_img_bin);


%% Real-Time simulation
% for i = 1:Nt
%     tic;
%     % 0: Load multi-echo volumes for current iteration
%     
%     for e = 1:Ne
%         fdyn_fn{e} = [wrf_me_fn{e} ',' num2str(i)]; % filename of dynamic functional image
%         currentVol{e} = spm_vol(fdyn_fn{e});
%         F_dyn_img{e} = spm_read_vols(currentVol{e}); % this is the unprocessed image
%         F_dyn{e}(:,i) = F_dyn_img{e}(:);
%         F_dyn_resliced{e}(:,i) = F_dyn{e}(:,i);
%         F_dyn_resliced_masked{e}(I_mask,i) = F_dyn_resliced{e}(I_mask,i);
%     end
%     
%     
%     % 2: estimate T2* and S0 per volume
%     S_pv = [F_dyn_resliced{1}(I_mask,i)'; F_dyn_resliced{2}(I_mask,i)'; F_dyn_resliced{3}(I_mask,i)'];
%     S_pv = max(S_pv,1e-11); %negative or zero signal values not possible
%     b_pv = X\log(S_pv);
%     S0_pv(I_mask,i)=exp(b_pv(1,:));
%     T2star_pv(I_mask,i)=1./b_pv(2,:);
%     
%     if isnan(b_pv)
%         disp(['isnan: i = ' i ])
%     end
%     
%     T2star_pv_corrected(I_mask,i) = T2star_pv(I_mask,i);
%     T2star_pv_corrected((T2star_pv_corrected(:,i)<0), i) = 0;
%     T2star_pv_corrected((T2star_pv_corrected(:,i)>100), i) = 0;
%     
%     
%     % 3: real-time combination using various methods
%     % combined_img2D = combine_echoes_v2(TE, method, weight_img, F, I_mask, N);
%     F = {F_dyn_resliced{1}(:,i), F_dyn_resliced{2}(:,i), F_dyn_resliced{3}(:,i)};
%     %     disp('before combine 1')
%     weight_img = T2star_pre_img;
%     combined_t2s_pre(:, i) = combine_echoes_v2(TE, 1, weight_img, F, I_mask, N);
%     %     disp('before combine 2')
%     weight_img = F_tSNR;
%     combined_tsnr_pre(:, i) = combine_echoes_v2(TE, 2, weight_img, F, I_mask, N);
%     %     disp('before combine 3')
%     weight_img = reshape(T2star_pv_corrected(:,i), Ni, Nj, Nk);
%     combined_t2s_rt(:, i) = combine_echoes_v2(TE, 1, weight_img, F, I_mask, N);
%     
%     %     % 2: GLM denoising
%     %     for e = 1:Ne
%     %         x_design = X_design(1:i, :);
%     %         beta_func = x_design\F_dyn_smoothed{e}(I_mask,1:i)'; % func = X*beta + e ==> beta = X\func ==> func_detrended = mp - X(i)*beta(i)
%     %         F_denoised{e} = F_dyn_smoothed{e}(I_mask,1:i)' - x_design(:, 2:(end-1))*beta_func(2:(end-1), :); % remove effects of all regressors except constant
%     %         F_denoised{e} = F_denoised{e}';
%     %         F_dyn_denoised{e}(I_mask,i) = F_denoised{e}(:,i);
%     % %         output = detrend4D(rf_me_fn{echo});
%     % %         S_preestimation_denoised{echo} = output.F_2D_detrended;
%     %     end
%     
%     T(i) = toc;
%     disp(['i=' num2str(i) ': ' num2str(T(i))]);
%     
% end

%% cluster analysis

tmap_fn = [stats_dir filesep 'spmT_0001.nii'];
xSPM = load([stats_dir filesep 'xSPM.mat']);
threshold = xSPM.xSPM.u;
[clusters, num] = findClustersMNI(tmap_fn, threshold);
Tind = 0;
max_val = 0;
val = 0;
for c = 1:num
    if clusters{c,2} > val
        val = clusters{c,2};
        Tind = c;
    end
end

I_cluster = clusters{Tind,1}(:,1);
tmap = spm_read_vols(spm_vol(tmap_fn));
tmap_reshaped = tmap(:);

cluster_map = zeros(Ni,Nj,Nk);
cluster_map_bin = cluster_map;
cluster_map = cluster_map(:);
cluster_map_bin = cluster_map_bin(:);
cluster_map(I_cluster) = tmap_reshaped(I_cluster);
cluster_map_bin(I_cluster) = 1;
cluster_map_img = reshape(cluster_map, Ni,Nj,Nk);
cluster_map_bin_img = reshape(cluster_map_bin, Ni,Nj,Nk);
cluster_threshold = val+2;
cluster_colormap = 'hot';
montage1 = createMontage(cluster_map_img, 8, 1, 'High-T-value cluster', cluster_colormap, 0, [0 cluster_threshold]);
%%
f = displayMaskContour(tmap, cluster_map_bin_img, 0, 5);

%%
srf = '/Users/jheunis/Desktop/rtme_tests/srsub-01_task-fingertapping.nii';
F = spm_read_vols(spm_vol(srf));
F_reshaped = reshape(F, Ni*Nj*Nk, Ntt);
clusterTL1 = F_reshaped(I_cluster,:);
clusterTL1ave = mean(clusterTL1,1);
clusterTL1ave_detrended = detrend(clusterTL1ave);
%%
figure; plot(1:Ntt, clusterTL1ave_detrended)
hold on; plot(1:Ntt, 10*convolved_task_design)

%%

% %%
% fn1 = '/Users/jheunis/Desktop/rtme_tests/rBroca_44.nii';
% fn2 = '/Users/jheunis/Desktop/rtme_tests/rBroca_45.nii';
%
% fn1_img = spm_read_vols(spm_vol(fn1));
% fn2_img = spm_read_vols(spm_vol(fn2));
% struct_img = spm_read_vols(spm_vol(preproc_data.wrstructural_fn));
%
% fn1_imgleft = fn1_img;
% fn1_imgleft(40:end, :, :) = 0;
% fn2_imgleft = fn2_img;
% fn2_imgleft(40:end, :, :) = 0;
%
% broca_mask1 = zeros(size(fn1_img));
% broca_mask1 = broca_mask1(:);
% broca_mask2 = broca_mask1;
%
% broca_mask1(find(fn1_imgleft(:)>0.2)) = 1;
% broca_mask2(find(fn2_imgleft(:)>0.2)) = 1;
%
% broca_mask = broca_mask1 | broca_mask2;
% I_broca = find(broca_mask);
% I_cluster = find(broca_mask);
%
% broca_mask_img = reshape(broca_mask, size(fn1_img));
%
% montage1 = createMontage(fn1_img, 8, 1, 'Broca 1', cluster_colormap, 0, [0 1]);
% montage2 = createMontage(fn2_img, 8, 1, 'Broca 2', cluster_colormap, 0, [0 1]);
% montage4 = createMontage(broca_mask_img, 8, 1, 'Broca 2', cluster_colormap, 0, [0 1]);
% montage3 = createMontage(struct_img, 8, 1, 'Struct', cluster_colormap, 0, [0 3500]);

%%






do_figs = 0;
if do_figs
    sub1_tSNR_ech02_mean = mean(sub1_tSNR_ech02 ,2);
    sub1_tSNR_preT2starcomb_mean = mean(sub1_tSNR_preT2starcomb ,2);
    sub1_tSNR_preTSNRcomb_mean = mean(sub1_tSNR_preTSNRcomb ,2);
    sub1_tSNR_rtT2starcomb_mean = mean(sub1_tSNR_rtT2starcomb ,2);
    sub1_perc_preT2starcomb_mean = mean(sub1_perc_preT2starcomb ,2);
    sub1_perc_preTSNRcomb_mean = mean(sub1_perc_preTSNRcomb ,2);
    sub1_perc_rtT2starcomb_mean = mean(sub1_perc_rtT2starcomb ,2);
    
    sub1_tSNR_preT2starcomb_mean(sub1_tSNR_preT2starcomb_mean>=250) = -100;
    sub1_tSNR_rtT2starcomb_mean(sub1_tSNR_rtT2starcomb_mean>=250) = -100;
    
    tsnr_threshold = 280;
    perc_threshold = 200;
    tsnr_colormap = 'hot';
    perc_colormap = 'parula';
    montage1 = createMontage(reshape(sub1_tSNR_ech02_mean, Ni, Nj, Nk), 8, 1, 'Ave tSNR - echo2', tsnr_colormap, 0, [0 tsnr_threshold]);
    montage2 = createMontage(reshape(sub1_tSNR_preT2starcomb_mean, Ni, Nj, Nk), 8, 1, 'Ave tSNR - preT2star comb', tsnr_colormap, 0, [0 tsnr_threshold]);
    montage3 = createMontage(reshape(sub1_tSNR_preTSNRcomb_mean, Ni, Nj, Nk), 8, 1, 'Ave tSNR - preTSNR comb', tsnr_colormap, 0, [0 tsnr_threshold]);
    montage4 = createMontage(reshape(sub1_tSNR_rtT2starcomb_mean, Ni, Nj, Nk), 8, 1, 'Ave tSNR - rtT2star comb', tsnr_colormap, 0, [0 tsnr_threshold]);
    montage5 = createMontage(reshape(sub1_perc_preT2starcomb_mean, Ni, Nj, Nk), 8, 1, 'Ave perc increase - preT2star comb', perc_colormap, 0, [0 perc_threshold]);
    montage6 = createMontage(reshape(sub1_perc_preTSNRcomb_mean, Ni, Nj, Nk), 8, 1, 'Ave perc increase - preTSNR comb', perc_colormap, 0, [0 perc_threshold]);
    montage7 = createMontage(reshape(sub1_perc_rtT2starcomb_mean, Ni, Nj, Nk), 8, 1, 'Ave perc increase - rtT2star comb', perc_colormap, 0, [0 perc_threshold]);
    
    numbins = 50;
    numvox = 25000;
    figure;
    subplot(2,2,1); histogram( sub1_tSNR_ech02_mean, numbins); axis([10 250 0 numvox]);
    subplot(2,2,2);histogram( sub1_tSNR_preT2starcomb_mean, numbins);axis([10 250 0 numvox]);
    subplot(2,2,3);histogram( sub1_tSNR_preTSNRcomb_mean, numbins);axis([10 250 0 numvox]);
    subplot(2,2,4);histogram( sub1_tSNR_rtT2starcomb_mean, numbins);axis([10 250 0 numvox]);
    figure;
    h1 = histogram(sub1_tSNR_ech02_mean, 'NumBins', numbins, 'FaceColor','r','EdgeColor','k','facealpha',0.5,'edgealpha',0.5);
    hold on;
    h2 = histogram(sub1_tSNR_preT2starcomb_mean, 'NumBins', numbins, 'FaceColor','g','EdgeColor','k','facealpha',0.5,'edgealpha',0.5);
    h3 = histogram( sub1_tSNR_preTSNRcomb_mean, 'NumBins', numbins, 'FaceColor','b','EdgeColor','k','facealpha',0.5,'edgealpha',0.5);
    h4 = histogram( sub1_tSNR_rtT2starcomb_mean, 'NumBins', numbins, 'FaceColor','c','EdgeColor','k','facealpha',0.5,'edgealpha',0.5);
    axis([10 250 0 numvox]);
    legend('Echo 2', 'Combined - preT2star', 'Combined - preTSNR', 'Combined - rtT2star');
end
%%









%
% kaas = sub1_tSNR_preT2starcomb_mean;
% kaas(kaas>=250) = -100;
% figure; hist( kaas );axis([0 250 0 1000000])
% figure; plot( kaas, '.k')
%






% %% Estimate parameters from averaged denoised dataset, or average estimated parameters
%
% % estimate from mean of denoised time series
% S0_mean = zeros(Nvox,1);
% T2star_mean = zeros(Nvox,1);
% T2star_mean_corrected = T2star_mean;
%
% S_mean = [mean(F_denoised{1},2)'; mean(F_denoised{2},2)'; mean(F_denoised{3},2)'; mean(F_denoised{4},2)'];
% S_mean = max(S_mean,1e-11); %negative or zero signal values not possible
% b_mean = X\log(S_mean);
% S0_mean(I_mask,:)=exp(b_mean(1,:));
% T2star_mean(I_mask,:)=1./b_mean(2,:);
%
% T2star_mean_corrected(I_mask) = T2star_mean(I_mask);
% T2star_mean_corrected((T2star_mean_corrected(:)<0)) = 0;
% T2star_mean_corrected((T2star_mean_corrected(:)>500)) = 0;
%
% T2star_mean_img = reshape(T2star_mean_corrected, Ni, Nj, Nk);
% S0_mean_img = reshape(S0_mean, Ni, Nj, Nk);
%
% % estimate per volume (pv), then average
% S0_pv = zeros(Nvox,Nt);
% T2star_pv = zeros(Nvox,Nt);
% T2star_pv_corrected = T2star_pv;
%
% for i = 1:Nt
% S_pv = [F_denoised{1}(:,i)'; F_denoised{2}(:,i)'; F_denoised{3}(:,i)'; F_denoised{1}(:,i)'];
% S_pv = max(S_pv,1e-11); %negative or zero signal values not possible
% b_pv = X\log(S_pv);
% S0_pv(I_mask,i)=exp(b_pv(1,:));
% T2star_pv(I_mask,i)=1./b_pv(2,:);
%
% if isnan(b_pv)
%     disp(['isnan: i = ' i ])
% end
% T2star_pv_corrected(I_mask,i) = T2star_pv(I_mask,i);
% T2star_pv_corrected((T2star_pv_corrected(:,i)<0), i) = 0;
% T2star_pv_corrected((T2star_pv_corrected(:,i)>500), i) = 0;
% end
%
% T2star_pv_mean = mean(T2star_pv_corrected, 2);
% S0_pv_mean = mean(S0_pv, 2);
%
% T2star_pv_mean_img = reshape(T2star_pv_mean, Ni, Nj, Nk);
% S0_pv_mean_img = reshape(S0_pv_mean, Ni, Nj, Nk);
% %%
%
% montage_m = createMontage(T2star_mean_img, 5, 1, 'T2star (estimated from temporal mean of timeseries)', 'hot', 0, [0 50]);
% montage_pvm = createMontage(T2star_pv_mean_img, 5, 1, 'T2star (estimated per volume, then averaged over time)', 'hot', 0, [0 50]);
%
%
% %% Figure
%
% figure;
% slice = 10;
% ax1 = subplot(2,4,1); imagesc(squeeze(T2star_img(:,:,slice,1)));
% ax2 = subplot(2,4,2); imagesc(squeeze(T2star_img(:,:,slice,1)));
% ax3 = subplot(2,4,3); imagesc(squeeze(T2star_img(:,:,slice,1)));
% ax4 = subplot(2,4,4); imagesc(squeeze(T2star_img(:,:,slice,1)));
% ax5 = subplot(2,4,5); imagesc(squeeze(T2star_img(:,:,slice,1)));
% ax6 = subplot(2,4,6); imagesc(squeeze(T2star_img(:,:,slice,1)));
% ax7 = subplot(2,4,7); imagesc(squeeze(T2star_img(:,:,slice,1)));
% ax8 = subplot(2,4,8); imagesc(squeeze(T2star_img(:,:,slice,1)));
%
% counter = zeros(Nt,10);
% % T2star_PREmean = mean(T2star_corrected(:, 150:end), 2);
% % T2star_PREmean_img = reshape(T2star_PREmean, Ni, Nj, Nk);
% % save('T2star_PREmean_img.mat', 'T2star_PREmean_img')
% % load('T2star_PREmean_img')
% %% real-time calcs and visuals using specific T2star map
%
% % T2star_PREmean = mean(T2star_corrected(:, 150:end), 2);
% % T2star_PREmean_img = reshape(T2star_PREmean, Ni, Nj, Nk);
% % save('T2star_PREmean_img.mat', 'T2star_PREmean_img')
% % load('T2star_PREmean_img')
% for i = 1:Nt
%     tic;
%
%     T2star_pv_img(:,:,:,i) = reshape(T2star_pv_corrected(:,i), Ni, Nj, Nk);
%     S0_pv_img(:,:,:,i) = reshape(S0(:,i), Ni, Nj, Nk);
%
%     T1_img = reshape(F_dyn_resliced_masked{1}(:,i), Ni, Nj, Nk);
%     T2_img = reshape(F_dyn_resliced_masked{2}(:,i), Ni, Nj, Nk);
%     T3_img = reshape(F_dyn_resliced_masked{3}(:,i), Ni, Nj, Nk);
%     T4_img = reshape(F_dyn_resliced_masked{4}(:,i), Ni, Nj, Nk);
%
%     T5_img = reshape(F_dyn_denoised{1}(:,i), Ni, Nj, Nk);
%     T6_img = reshape(F_dyn_denoised{2}(:,i), Ni, Nj, Nk);
%     T7_img = reshape(F_dyn_denoised{3}(:,i), Ni, Nj, Nk);
%     T8_img = reshape(F_dyn_denoised{4}(:,i), Ni, Nj, Nk);
%
%     S_combined_img(:,:,:,i) = combine_echoes(TE, T2star_mean_img, I_mask, 1, T1_img, T2_img, T3_img, T4_img);
%     S_combined_img2(:,:,:,i) = combine_echoes(TE, T2star_mean_img, I_mask, 1, T5_img, T6_img, T7_img, T8_img);
% %     S_combined_img(:,:,:,i) = combine_echoes(TE, T2star_PREmean_img, I_mask, 1, T1_img, T2_img, T3_img);
% %     S_combined_img2(:,:,:,i) = combine_echoes(TE, T2star_PREmean_img, I_mask, 1, T4_img, T5_img, T6_img);
%
%
%     imagesc(ax1, squeeze(T1_img(:,:,slice))); colormap(ax1, 'bone'); colorbar(ax1); %ax1.CLim = [0 15e5];
%     imagesc(ax2, squeeze(T2_img(:,:,slice))); colormap(ax2, 'bone'); colorbar(ax2); %ax2.CLim = [0 15e5];
%     imagesc(ax3, squeeze(T3_img(:,:,slice))); colormap(ax3, 'bone'); colorbar(ax3); %ax3.CLim = [0 15e5];
%     imagesc(ax4, squeeze(T4_img(:,:,slice))); colormap(ax4, 'bone'); colorbar(ax4); %ax4.CLim = [0 15e5];
%
%     imagesc(ax5, squeeze(S_combined_img(:,:,slice,i))); colormap(ax5, 'bone'); colorbar(ax5); %ax5.CLim = [0 15e5];
%     imagesc(ax6, squeeze(S_combined_img2(:,:,slice,i))); colormap(ax6, 'bone'); colorbar(ax6); %ax6.CLim = [0 15e5];
%
%     imagesc(ax7, squeeze(T2star_pv_img(:,:,slice,i))); colormap(ax7, 'gray'); colorbar(ax7); %ax7.CLim = [0 300];
%     imagesc(ax8, squeeze(S0_pv_img(:,:,slice,i))); colormap(ax7, 'gray'); colorbar(ax8); %ax8.CLim = [0 15e5];
%
%     drawnow;
%     T(i) = toc;
%     disp(['i=' num2str(i) ': ' num2str(T(i))]);
% %     % Estimate maps
% %     for i = 1:Nmaskvox
% %         for j = 1:Nt
% %             % TODO: incorporate both smoothed and unsmoothed preprocessed
% %             % signal here. Or either.
% %             S1=[S_preestimation_denoised{1,1}(I_mask(i),j); S_preestimation_denoised{2,1}(I_mask(i),j); S_preestimation_denoised{3,1}(I_mask(i),j)];
% %             S1=max(S1,1e-11); %negative or zero signal values not possible
% %             b1=X\log(S1);
% %             S0(I_mask(i),j)=exp(b1(1));
% %             T2star(I_mask(i),j)=1/b1(2);
% %
% %             if isnan(b1)
% %                 disp(['isnan: i = ' i '; j = ' j])
% %             end
% %         end
% %     end
%
% end
%
% %% generate and plot tSNR images
% F_tsnrR = cell(Ne+3,1);
% F_tsnrS = cell(Ne+3,1);
% F_tsnrR_img = cell(Ne+3,1);
% F_tsnrS_img = cell(Ne+3,1);
%
% for e = 1:Ne
%     clear m1 stddev1 m2 stddev2;
%
%     m1 = mean(F_dyn_resliced_masked{e}(I_mask,:), 2);
%     stddev1 = std(F_dyn_resliced_masked{e}(I_mask,:), 0, 2);
%
%     m2 = mean(F_dyn_smoothed_masked{e}(I_mask,:), 2);
%     stddev2 = std(F_dyn_smoothed_masked{e}(I_mask,:), 0, 2);
%
%     F_tsnrR{e} = zeros(Nvox,1);
%     F_tsnrR{e}(I_mask) = m1./stddev1;
%     F_tsnrR{e}(isnan(F_tsnrR{e}))=0;
%     F_tsnrR_img{e} = reshape(F_tsnrR{e}, Ni, Nj, Nk);
%
%     F_tsnrS{e} = zeros(Nvox,1);
%     F_tsnrS{e}(I_mask) = m2./stddev2;
%     F_tsnrS{e}(isnan(F_tsnrS{e}))=0;
%     F_tsnrS_img{e} = reshape(F_tsnrS{e}, Ni, Nj, Nk);
%
% end
%
%
% combined = reshape(S_combined_img, Nvox, Nt);
% m = mean(combined(I_mask,:), 2);
% stddev = std(combined(I_mask,:), 0, 2);
% F_tsnrR{Ne+1} = zeros(Nvox,1);
% F_tsnrR{Ne+1}(I_mask) = m./stddev;
% F_tsnrR{Ne+1}(isnan(F_tsnrR{Ne+1}))=0;
% F_tsnrR_img{Ne+1} = reshape(F_tsnrR{Ne+1}, Ni, Nj, Nk);
%
% combined = reshape(S_combined_img2, Nvox, Nt);
% m = mean(combined(I_mask,:), 2);
% stddev = std(combined(I_mask,:), 0, 2);
% F_tsnrS{Ne+1} = zeros(Nvox,1);
% F_tsnrS{Ne+1}(I_mask) = m./stddev;
% F_tsnrS{Ne+1}(isnan(F_tsnrS{Ne+1}))=0;
% F_tsnrS_img{Ne+1} = reshape(F_tsnrS{Ne+1}, Ni, Nj, Nk);
%
% %% plot tSNR images
% figure;
% slice = 10;
% ax10 = subplot(3,2,1); imagesc(ax10, squeeze(F_tsnrR_img{1}(:,:,slice))); colormap(ax10, 'hot'); colorbar(ax10); ax10.CLim = [0 400];
% ax11 = subplot(3,2,2); imagesc(ax11, squeeze(F_tsnrR_img{2}(:,:,slice))); colormap(ax11, 'hot'); colorbar(ax11); ax11.CLim = [0 400];
% ax12 = subplot(3,2,3); imagesc(ax12, squeeze(F_tsnrR_img{3}(:,:,slice))); colormap(ax12, 'hot'); colorbar(ax12); ax12.CLim = [0 400];
% ax13 = subplot(3,2,4); imagesc(ax13, squeeze(F_tsnrR_img{4}(:,:,slice))); colormap(ax13, 'hot'); colorbar(ax13); ax13.CLim = [0 400];
% ax14 = subplot(3,2,5); imagesc(ax14, squeeze(F_tsnrR_img{5}(:,:,slice))); colormap(ax14, 'hot'); colorbar(ax14); ax14.CLim = [0 400];
%
% %%
% figure;
% slice = 10;
% ax10 = subplot(3,2,1); imagesc(ax10, squeeze(F_tsnrS_img{1}(:,:,slice))); colormap(ax10, 'hot'); colorbar(ax10); ax10.CLim = [0 400];
% ax11 = subplot(3,2,2); imagesc(ax11, squeeze(F_tsnrS_img{2}(:,:,slice))); colormap(ax11, 'hot'); colorbar(ax11); ax11.CLim = [0 400];
% ax12 = subplot(3,2,3); imagesc(ax12, squeeze(F_tsnrS_img{3}(:,:,slice))); colormap(ax12, 'hot'); colorbar(ax12); ax12.CLim = [0 400];
% ax13 = subplot(3,2,4); imagesc(ax13, squeeze(F_tsnrS_img{4}(:,:,slice))); colormap(ax13, 'hot'); colorbar(ax13); ax13.CLim = [0 400];
% ax14 = subplot(3,2,5); imagesc(ax14, squeeze(F_tsnrS_img{5}(:,:,slice))); colormap(ax14, 'hot'); colorbar(ax14); ax14.CLim = [0 400];
%
% %%
% clear montage1 montage2 montage3 montage4;
% montage1 = createMontage(F_tsnrR_img{2}, 4, 1, 'tSNR (realigned echo 2)', 'hot', 0, [0 250]);
% % montage2 = createMontage(F_tsnrS_img{2}, 4, 1, 'tSNR (smoothed echo 2)', 'hot', 0, [0 400]);
% montage3 = createMontage(F_tsnrR_img{5}, 4, 1, 'tSNR (realigned combined)', 'hot', 0, [0 250]);
% % montage4 = createMontage(F_tsnrS_img{5}, 4, 1, 'tSNR (smoothed combined)', 'hot', 0, [0 400]);
%
%
% %%
% % %%
% % % image = S_preestimation_denoised{2};
% % % S_combined_2D = reshape(S_combined_img, Ni*Nj*Nk, Nt);
% % % image = S_combined_2D;
% % T2star_2D = reshape(T2star_img, Ni*Nj*Nk, Nt);
% % image = T2star_2D;
% %
% % % Statistical measures
% % F2D_mean = mean(image, 2);
% % F2D_stddev = std(image, 0, 2);
% % F2D_zstat = (image - F2D_mean)./F2D_stddev;
% % F2D_zstat(isnan(F2D_zstat))=0;
% % F2D_zstat_mean = mean(abs(F2D_zstat),1);
% % Zstat_mean = mean(F2D_zstat_mean);
% % F2D_psc = 100*(image./repmat(F2D_mean, 1, Nt)) - 100;
% % F2D_psc(isnan(F2D_psc))=0;
% % F2D_diff = [zeros(1, Ni*Nj*Nk); diff(image')]';
% % F2D_DVARS = var(F2D_diff);
% % % tSNR
% % tSNR_2D = F2D_mean./F2D_stddev;
% % tSNR_2D(isnan(tSNR_2D))=0;
% % tSNR_brain = mean(tSNR_2D(I_mask));
% % tSNR_GM = mean(tSNR_2D(I_GM));
% % tSNR_WM = mean(tSNR_2D(I_WM));
% % tSNR_CSF = mean(tSNR_2D(I_CSF));
% % % Metrics
% % disp(['Mean Zscore: ' num2str(Zstat_mean)])
% % disp(['tSNR (brain): ' num2str(tSNR_brain)])
% % disp(['tSNR (GM): ' num2str(tSNR_GM)])
% % disp(['tSNR (WM): ' num2str(tSNR_WM)])
% % disp(['tSNR (CSF): ' num2str(tSNR_CSF)])
% %
% % % 3D and 4D images
% % mask_3D = reshape(mask_reshaped, Ni, Nj, Nk);
% % tSNR_3D = reshape(tSNR_2D, Ni, Nj, Nk);
% % F3D_mean = reshape(F2D_mean, Ni, Nj, Nk);
% % F3D_stddev = reshape(F2D_stddev, Ni, Nj, Nk);
% %
% % montage2 = createMontage(F3D_mean, 5, 1, 'Mean EPI (whole image)', 'gray');
% % montage3 = createMontage(F3D_stddev, 5, 1, 'Standard deviation (whole image)', 'parula');
% % montage1 = createMontage(tSNR_3D, 5, 1, 'tSNR (whole image)', 'hot');
%
% %% PSC
% % echo2_resliced = reshape(F_dyn_resliced_masked{2}, Ni, Nj, Nk, Nt);
% % echo2_denoised = reshape(F_dyn_denoised{2}, Ni, Nj, Nk, Nt);
% echo2_resliced = F_dyn_smoothed_masked{2};
% echo2_denoised = F_dyn_denoised{2};
% PSC_echo2r = getPercentageSignalChange(echo2_resliced, 6, 7, 16, [17;49;81;113;145;177], 208);
% PSC_echo2d = getPercentageSignalChange(echo2_denoised, 6, 7, 16, [17;49;81;113;145;177], 208);
%
% S1 = reshape(S_combined_img, Nvox, Nt);
% S2 = reshape(S_combined_img, Nvox, Nt);
% PSC_Sr = getPercentageSignalChange(S1, 6, 7, 16, [17;49;81;113;145;177], 208);
% PSC_Sd = getPercentageSignalChange(S2, 6, 7, 16, [17;49;81;113;145;177], 208);
%
% PSC_blocks_maskedr = zeros(Ni*Nj*Nk, Nt);
% PSC_blocks_maskedr(I_mask, :) = PSC_echo2r.PSC_blocks(I_mask, :);
% PSC_blocks_4Dr = reshape(PSC_blocks_maskedr, Ni, Nj, Nk, Nt);
%
% PSC_blocks_maskedd = zeros(Ni*Nj*Nk, Nt);
% PSC_blocks_maskedd(I_mask, :) = PSC_echo2d.PSC_blocks(I_mask, :);
% PSC_blocks_4Dd = reshape(PSC_blocks_maskedd, Ni, Nj, Nk, Nt);
%
% PSC_blocks_maskedSr = zeros(Ni*Nj*Nk, Nt);
% PSC_blocks_maskedSr(I_mask, :) = PSC_Sr.PSC_blocks(I_mask, :);
% PSC_blocks_4DSr = reshape(PSC_blocks_maskedSr, Ni, Nj, Nk, Nt);
%
% PSC_blocks_maskedSd = zeros(Ni*Nj*Nk, Nt);
% PSC_blocks_maskedSd(I_mask, :) = PSC_Sd.PSC_blocks(I_mask, :);
% PSC_blocks_4DSd = reshape(PSC_blocks_maskedSd, Ni, Nj, Nk, Nt);
%
% temporal_mask = ((PSC_echo2d.task_blocks==5)|(PSC_echo2d.task_blocks==6));
%
% PSC_ave_resliced = mean(PSC_blocks_maskedr(:,temporal_mask), 2);
% PSC_ave_resliced3D = reshape(PSC_ave_resliced, Ni, Nj, Nk);
% PSC_ave_denoised = mean(PSC_blocks_maskedd(:,temporal_mask), 2);
% PSC_ave_denoised3D = reshape(PSC_ave_denoised, Ni, Nj, Nk);
%
% PSC_Save_resliced = mean(PSC_blocks_maskedSr(:,temporal_mask), 2);
% PSC_Save_resliced3D = reshape(PSC_Save_resliced, Ni, Nj, Nk);
% PSC_Save_denoised = mean(PSC_blocks_maskedSd(:,temporal_mask), 2);
% PSC_Save_denoised3D = reshape(PSC_Save_denoised, Ni, Nj, Nk);
%
% figure;
% slice = 10;
% ax21 = subplot(2,2,1); imagesc(ax21, squeeze(PSC_ave_resliced3D(:,:,slice))); colormap(ax21, 'hot'); colorbar(ax21); ax21.CLim = [0 10];
% ax22 = subplot(2,2,2); imagesc(ax22, squeeze(PSC_ave_denoised3D(:,:,slice))); colormap(ax22, 'hot'); colorbar(ax22); ax22.CLim = [0 10];
% ax23 = subplot(2,2,3); imagesc(ax23, squeeze(PSC_Save_resliced3D(:,:,slice))); colormap(ax23, 'hot'); colorbar(ax23); ax23.CLim = [0 10];
% ax24 = subplot(2,2,4); imagesc(ax24, squeeze(PSC_Save_denoised3D(:,:,slice))); colormap(ax24, 'hot'); colorbar(ax24); ax24.CLim = [0 10];
%
%
%
% %%
% montage3 = createMontage(PSC_ave_denoised3D, 4, 1, 'PSC (denoised echo 2)', 'hot', 0, [-5 5]);
% montage4 = createMontage(PSC_Save_denoised3D, 4, 1, 'PSC (denoised combined)', 'hot', 0, [-5 5]);
%
%
%
%
% % PSC_all_masked = zeros(Ni*Nj*Nk, Nt);
% % PSC_all_masked(I_mask, :) = F2D_psc(I_mask, :);
% % PSC_all_4D = reshape(PSC_all_masked, Ni, Nj, Nk, Nt);
%
% % new_nii = make_nii(PSC_all_4D, [3.5 3.5 4.5], [27 37 17]);
% % save_nii(new_nii, 'T2S_PSC_all.nii')
% % new_nii = make_nii(PSC_blocks_4D, [3.5 3.5 4.5], [27 37 17]);
% % save_nii(new_nii, 'T2S_PSC_blocks.nii')
% %
% % draw_brain_views_gui([20 20 20], PSC_blocks_4D)
% % draw_brain_views_gui([20 20 20], PSC_all_4D)
%
%
%
% %% cluster analysis
% stats_dir = [data_dir filesep subj_dir filesep 'stats'];
% tmap_fn = [stats_dir filesep 'spmT_0001.nii'];
% xSPM = load([stats_dir filesep 'xSPM.mat']);
% threshold = xSPM.xSPM.u;
% SPM = load([stats_dir filesep 'SPM.mat']);
% % convolved_task_timecourse = SPM.SPM.xX.X(:,1);
%
% [clusters, num] = findClustersMNI(tmap_fn, threshold);
%
% % [Tmax, Tind] = max(clusters{:,2})
% Tind = 0;
% max_val = 0;
% val = 0;
% for c = 1:num
%     if clusters{c,2} > val
%         val = clusters{c,2};
%         Tind = c;
%     end
% end
%
% I_cluster = clusters{Tind,1}(:,1);
% tmap = spm_read_vols(spm_vol(tmap_fn));
% [Ni,Nj,Nk] = size(tmap);
% cluster_map = zeros(Ni,Nj,Nk);
% cluster_map = cluster_map(:);
% cluster_map(I_cluster) = 1;
% cluster_map_img = reshape(cluster_map, Ni,Nj,Nk);
% % cluster_nifti = [data_dir_new filesep subj_dir filesep 'main_cluster.nii'];
% % % disp('Creating cluster nifti  ...')
% % if ~exist(cluster_nifti, 'file')
% %     new_nii = make_nii(cluster_map_img, [3.5 3.5 4.5], [27 37 17]);
% %     save_nii(new_nii, cluster_nifti)
% % end
%
% %% display clusters
% fff = displayMaskContour(tmap, cluster_map_img, 0, 4);
%
%
% %%
% % print(montage1.f, 'tsnr_whole', '-dpng')
% % montage4 = createMontage(tSNR_3D_masked, 5, 1, 'tSNR (brain)', 'hot');
%
%
%
%
% % % The plot (with FD, DVARS, and mean Zscore per volume)
% % GM_img = F2D_psc(I_GM, :);
% % WM_img = F2D_psc(I_WM, :);
% % CSF_img = F2D_psc(I_CSF, :);
% % all_img = [GM_img; WM_img; CSF_img];
% % line1_pos = numel(I_GM);
% % line2_pos = numel(I_GM) + numel(I_WM);
% % tf = figure;
% % fontsizeL = 14;
% % fontsizeM = 11;
% % ax1 = subplot(7,1,4:7);
% % imagesc(ax1, all_img); colormap(gray); caxis(intensity_scale);
% % title(ax1, 'thePlotSpm','fontsize',fontsizeL)
% % ylabel(ax1, 'Voxels','fontsize',fontsizeM)
% % xlabel(ax1, 'fMRI volumes','fontsize',fontsizeM)
% % hold on; line([1 Nt],[line1_pos line1_pos],  'Color', 'b', 'LineWidth', 2 )
% % line([1 Nt],[line2_pos line2_pos],  'Color', 'r', 'LineWidth', 2 )
% % hold off;
% % ax2 = subplot(7,1,1);
% % plot(ax2, FD_measures.FD, 'LineWidth', 2); grid;
% % set(ax2,'Xticklabel',[]);
% % title(ax2, 'FD','fontsize',fontsizeL)
% % ylabel(ax2, 'mm','fontsize',fontsizeM)
% % ax3 = subplot(7,1,2);
% % plot(ax3, DVARS, 'LineWidth', 2); grid;
% % set(ax3,'Xticklabel',[]);
% % title(ax3, 'DVARS','fontsize',fontsizeL)
% % ylabel(ax3, 'a.u.','fontsize',fontsizeM)
% % ax4 = subplot(7,1,3);
% % plot(ax4, F2D_zstat_mean, 'LineWidth', 2); grid;
% % set(ax4,'Xticklabel',[]);
% % title(ax4, 'Z-score','fontsize',fontsizeL)
% % ylabel(ax4, 'a.u.','fontsize',fontsizeM)
% % print(tf, 'timeseries_summary', '-dpng')
%
%
%
%
%
%
%
%
%
% %% PIPELINE 3: COMBINED ECHO PROCESSING (STANDARD GLM DENOISING OF COMBINED ECHO IMAGES)
% disp('PIPELINE 3: COMBINED ECHO PROCESSING ...')
%
% disp('Combining smoothed multi-echo images ...')
% S_combined_img = combine_echoes(TE, T2star_img, I_mask, 1, T1_img, T2_img, T3_img);
% % S_combined_img = combine_echoes([echo time vector], T2star_img, I_overlap, [method (1 = timeseries, 2 = average, 3 = average blocks)], T1_img, T2_img, T3_img);
%
% %Then create new nifti
% disp('Creating cleaned combined nifti  ...')
% if ~exist(analysis_nii_fn{3}, 'file')
%     new_nii = make_nii(S_combined_img, [3.5 3.5 4.5], [27 37 17]);
%     save_nii(new_nii, analysis_nii_fn{3})
% end
%
% %Then run stats
% disp('Running SPM stats on cleaned combined nifti ...')
% if ~exist([mni_stats_dir{3} filesep 'spmT_0001.nii'], 'file')
%     runSpmStats(jobs_dir, mni_stats_dir{3}, analysis_nii_fn{3}, Nt, Nregr, regressors_fn{2})
%     pause(5);
%     save([mni_stats_dir{3} filesep 'xSPM.mat'], 'xSPM')
% end
%
% %% correlation stuff
% disp('Correlation calcs')
% tmap_fn = ['/Users/jheunis/Documents/MATLAB/ME/Volunteer_data/V12all/' subj_dir '/Task/Pipeline 2 - TE2/MNIstats/spmT_0001.nii'];
% tmap = spm_read_vols(spm_vol(tmap_fn));
% load(['/Users/jheunis/Documents/MATLAB/ME/Volunteer_data/V12all/' subj_dir '/Task/Pipeline 2 - TE2/MNIstats/xSPM.mat'])
% tmap_thresholded = (tmap > xSPM.u);
% tmap_thresholded0 = (tmap > 0);
%
% [vals, inds] = sort(tmap(:), 'descend');
% vals_thresholded = vals > xSPM.u;
% Nu = numel(find(vals_thresholded));
% vals_thresholded0 = vals > 0;
% N0 = numel(find(vals_thresholded0));
%
% [r0, p0] = corr(T2star(inds(1:N0), :)', S0(inds(1:N0), :)');
%
% r_diag0 = zeros(1,N0);
% for j = 1:numel(r_diag0)
%     r_diag0(j) = r0(j,j);
% end
% new_vals = flipud(vals(1:N0));
% new_r = fliplr(r_diag0);
% figure; plot(new_vals, new_r, '.r'); grid; hold on;
% line([xSPM.u xSPM.u], [-1 1]);
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% title(['Subject-' num2str(subj) '-' num2str(run) ': Correlation between T2star and S0 (per voxel) for Tmap voxels where T>0 (TE2)']);
% xlabel(['T-values from Tmap (SPM threshold at ' num2str(xSPM.u) ')']);
% ylabel('Pearson coefficient R per voxel');
%
%
%
% disp('Correlation calcs 2')
% tmap_fn = ['/Users/jheunis/Documents/MATLAB/ME/Volunteer_data/V12all/' subj_dir '/Task/Pipeline 3 - Combined_TE/MNIstats/spmT_0001.nii'];
% tmap = spm_read_vols(spm_vol(tmap_fn));
% load(['/Users/jheunis/Documents/MATLAB/ME/Volunteer_data/V12all/' subj_dir '/Task/Pipeline 3 - Combined_TE/MNIstats/xSPM.mat'])
% tmap_thresholded = (tmap > xSPM.u);
% tmap_thresholded0 = (tmap > 0);
%
% [vals, inds] = sort(tmap(:), 'descend');
% vals_thresholded = vals > xSPM.u;
% Nu = numel(find(vals_thresholded));
% vals_thresholded0 = vals > 0;
% N0 = numel(find(vals_thresholded0));
%
% [r0, p0] = corr(T2star(inds(1:N0), :)', S0(inds(1:N0), :)');
%
% r_diag0 = zeros(1,N0);
% for j = 1:numel(r_diag0)
%     r_diag0(j) = r0(j,j);
% end
% new_vals = flipud(vals(1:N0));
% new_r = fliplr(r_diag0);
% figure; plot(new_vals, new_r, '.r'); grid; hold on;
% line([xSPM.u xSPM.u], [-1 1]);
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% title(['Subject-' num2str(subj) '-' num2str(run) ': Correlation between T2star and S0 (per voxel) for Tmap voxels where T>0 (Combined_TE)']);
% xlabel(['T-values from Tmap (SPM threshold at ' num2str(xSPM.u) ')']);
% ylabel('Pearson coefficient R per voxel');

