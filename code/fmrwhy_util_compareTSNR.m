function fmrwhy_util_compareTSNR(tsnr_fns, mask_fn, roi_fns, compare_roi_txt, tsnr_saveAs_fns, perc_saveAs_fns, distr_saveAs_fn)


% Settings
percentage_threshold = 160;
tsnr_threshold = 320;

% Get mask indices
[p_mask, frm, rg, dim_mask] = fmrwhy_util_readOrientNifti(mask_fn);
mask_img = p_mask.nii.img;
I_mask = find(mask_img(:));

% Grab tsnr image data and calculate percentage difference of all images in the array from first image in the array
tsnr_img = {};
diff = {};
perc_diff = {};
for i = 1:numel(tsnr_fns)
    [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(tsnr_fns{i});
%    tsnr_img{i} = double(p.nii.img);
    tsnr_img{i} = fmrwhy_util_maskImage(double(p.nii.img), mask_img);
    if i > 1
        diff{i-1} = tsnr_img{i} - tsnr_img{1};
        perc_diff{i-1} = diff{i-1}./tsnr_img{1}*100;
        perc_diff{i-1}(perc_diff{i-1}>percentage_threshold) = percentage_threshold;
    end
end

% Grab+load ROI image data; get ROI indices; combine ROI image data into a single overlay image
roi_img = {};
I_roi = {};
overlay_img = zeros(dim_mask);
for i = 1:numel(roi_fns)
    [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(roi_fns{i});
    roi_img{i} = fmrwhy_util_createBinaryImg(p.nii.img, 0.1);
    I_roi{i} = find(roi_img{i}(:));
    overlay_img = overlay_img | roi_img{i};
end

% ---
% Output 1: Montage plots of tsnr with ROI overlays
% ---
for i = 1:numel(tsnr_img)
    if ~exist(tsnr_saveAs_fns{i}, 'file')
        overlaymontage = fmrwhy_util_createOverlayMontage(tsnr_img{i}, overlay_img, 9, 1, '', 'hot', 'off', 'max', [0 250], [33, 168, 10], tsnr_saveAs_fns{i});
    end
end

% ---
% Output 2: Montage plots of perc diff with ROI overlays
% ---
for i = 1:numel(perc_diff)
    if ~exist(perc_saveAs_fns{i}, 'file')
        overlaymontage = fmrwhy_util_createOverlayMontage(perc_diff{i}, overlay_img, 9, 1, '', 'parula', 'off', 'max', [0 300], [215,25,28], perc_saveAs_fns{i});
    end
end

% ---
% Output 3: Raincloud plot
% ---
[d,f,e] = fileparts(distr_saveAs_fn);

[cb] = cbrewer('qual','Set3',12,'pchip');
cl(1,:) = cb(4,:);
cl(2,:) = cb(1,:);

f1 = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
% tSNR distribution in whole brain
ax1 = subplot(1,1,1);
h1= raincloud_plot('X', tsnr_img{1}(I_mask) , 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .2, 'dot_dodge_amount', .2, 'box_col_match', 0);
h2= raincloud_plot('X', tsnr_img{2}(I_mask), 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .4, 'dot_dodge_amount', .4, 'box_col_match', 0);
h3= raincloud_plot('X', tsnr_img{3}(I_mask), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .6, 'dot_dodge_amount', .6, 'box_col_match', 0);
h4= raincloud_plot('X', tsnr_img{4}(I_mask), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .8, 'dot_dodge_amount', .8, 'box_col_match', 0);
legend([h1{1} h2{1} h3{1} h4{1}], {'Echo 2', 'T2{\ast}-combined', 'tSNR-combined', 'TE-combined'})
title('tSNR distributions within whole brain mask')
set(gca,'XLim', [-10 tsnr_threshold], 'YLim', [-.02 .02]);
hold(ax1, 'on')
plot(ax1, [0 0],[0 -0.2],'k')
hold(ax1, 'off')
box off

brain_saveAs_fn = [d filesep f '_brain' e];
if ~exist(brain_saveAs_fn, 'file')
    print(f1, brain_saveAs_fn,'-dpng', '-r0');
end

f2 = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
% tSNR distribution in ROI 1
%ax2 = subplot(3,2,3)
ax2 = subplot(2,2,1);
h1= raincloud_plot('X', tsnr_img{1}(I_roi{1}) , 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
h2= raincloud_plot('X', tsnr_img{2}(I_roi{1}), 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
h3= raincloud_plot('X', tsnr_img{3}(I_roi{1}), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
h4= raincloud_plot('X', tsnr_img{4}(I_roi{1}), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0);
legend([h1{1} h2{1} h3{1} h4{1}], {'Echo 2', 'T2{\ast}-combined', 'tSNR-combined', 'TE-combined'})
title(['tSNR distributions within ' compare_roi_txt{1}] )
set(gca,'XLim', [-10 tsnr_threshold], 'YLim', [-.02 .02]);
hold(ax2, 'on')
plot(ax2, [0 0],[0 -0.2],'k')
hold(ax2, 'off')
box off

% tSNR distribution in ROI 2
%ax3 = subplot(3,2,5)
ax3 = subplot(2,2,3);
h1= raincloud_plot('X', tsnr_img{1}(I_roi{2}) , 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
h2= raincloud_plot('X', tsnr_img{2}(I_roi{2}), 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
h3= raincloud_plot('X', tsnr_img{3}(I_roi{2}), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
h4= raincloud_plot('X', tsnr_img{4}(I_roi{2}), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0);
legend([h1{1} h2{1} h3{1} h4{1}], {'Echo 2', 'T2{\ast}-combined', 'tSNR-combined', 'TE-combined'})
title(['tSNR distributions within ' compare_roi_txt{2}] )
set(gca,'XLim', [-10 tsnr_threshold], 'YLim', [-.02 .02]);
hold(ax3, 'on')
plot(ax3, [0 0],[0 -0.2],'k')
hold(ax3, 'off')
box off

% Distribution of percentage difference in tSNR in ROI 1
%ax4 = subplot(3,2,4)
ax4 = subplot(2,2,2);
h1= raincloud_plot('X', perc_diff{1}(I_roi{1}) , 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
h2= raincloud_plot('X', perc_diff{2}(I_roi{1}), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
h3= raincloud_plot('X', perc_diff{3}(I_roi{1}), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
legend([h1{1} h2{1} h3{1}], {'T2{\ast}-combined', 'tSNR-combined', 'TE-combined'})
title(['Percentage difference distributions within ' compare_roi_txt{1}] )
set(gca,'XLim', [-50 percentage_threshold], 'YLim', [-.02 .04]);
hold(ax4, 'on')
plot(ax4, [0 0],[0 -0.4],'k')
hold(ax4, 'off')
box off

% Distribution of percentage difference in tSNR in ROI 2
%ax5 = subplot(3,2,6)
ax5 = subplot(2,2,4);
h1= raincloud_plot('X', perc_diff{1}(I_roi{2}) , 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
h2= raincloud_plot('X', perc_diff{2}(I_roi{2}), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
h3= raincloud_plot('X', perc_diff{3}(I_roi{2}), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
legend([h1{1} h2{1} h3{1}], {'T2{\ast}-combined', 'tSNR-combined', 'TE-combined'})
title(['Percentage difference distributions within ' compare_roi_txt{2}] )
set(gca,'XLim', [-50 percentage_threshold], 'YLim', [-.02 .04]);
hold(ax5, 'on')
plot(ax5, [0 0],[0 -0.4],'k')
hold(ax5, 'off')
box off

roi_saveAs_fn = [d filesep f '_roi' e];
if ~exist(roi_saveAs_fn, 'file')
    print(f2, roi_saveAs_fn,'-dpng', '-r0');
end
