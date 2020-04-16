function fmrwhy_util_compareTSNR(tsnr_fns, mask_fn, roi_fns, I_roi_distr, roi_text, tsnr_saveAs_fns, perc_saveAs_fns, distr_saveAs_fn)

% I_roi_distr is for distribution plot of percdiff in task roi. send 0 if rest.
% TODO 16 april 2020: this needs to be updated such that I_roi_distr is not sent as a parameter
% but rather as the task roi filename, such that correct read/write functions can be used inside this function
% to ensure consistency or nifti orientation.


%[p1, frm1, rg1, dim1] = fmrwhy_util_readNifti(t2star_fn);
%[p2, frm2, rg2, dim2] = fmrwhy_util_readNifti(s0_fn);
%t2star_montage = fmrwhy_util_createMontage(p1.nii.img, 9, 1, 'T2star', 'hot', 'on', 'max');
%colorbar; % caxis([0 200]);
%s0_montage = fmrwhy_util_createMontage(p2.nii.img, 9, 1, 'S0', 'parula', 'on', 'max');
%colorbar;
percentage_threshold = 200;

% Get mask indices
[p_mask, frm, rg, dim_mask] = fmrwhy_util_readNifti(mask_fn);
mask_img = p_mask.nii.img;
I_mask = find(mask_img(:));

% Grab tsnr image data and calculate percentage difference from first tsnr image
tsnr_img = {};
diff = {};
perc_diff = {};
for i = 1:numel(tsnr_fns)
    [p, frm, rg, dim] = fmrwhy_util_readNifti(tsnr_fns{i});
    tsnr_img{i} = p.nii.img .* mask_img;
    if i > 1
        diff{i-1} = tsnr_img{i} - tsnr_img{1};
        perc_diff{i-1} = diff{i-1}./tsnr_img{1}*100;
        perc_diff{i-1}(perc_diff{i-1}>percentage_threshold) = percentage_threshold;
    end
end

% Grab roi image data; combined them if more than one
roi_img = zeros(dim_mask);
overlay_img = {};
for i = 1:numel(roi_fns)
    [p, frm, rg, dim] = fmrwhy_util_readNifti(roi_fns{i});
    overlay_img{i} = fmrwhy_util_createBinaryImg(p.nii.img, 0.1);
    roi_img = roi_img | overlay_img{i};
end
%I_roi = find(roi_img(:));
%
%assignin('base','overlay_img',overlay_img)
%assignin('base','perc_diff',perc_diff)
%assignin('base','tsnr_img',tsnr_img)
%assignin('base','I_mask',I_mask)


% Montage plots of tsnr and perc diff
for i = 1:numel(tsnr_img)

    if ~exist(tsnr_saveAs_fns{i}, 'file')
        overlaymontage = fmrwhy_util_createOverlayMontage(tsnr_img{i}, roi_img, 9, 1, '', 'hot', 'off', 'max', [0 200], tsnr_saveAs_fns{i});
    end
end
for i = 1:numel(perc_diff)
    if ~exist(perc_saveAs_fns{i}, 'file')
        overlaymontage = fmrwhy_util_createOverlayMontage(perc_diff{i}, roi_img, 9, 1, '', 'parula', 'off', 'max', [0 300], perc_saveAs_fns{i});
    end
end


%% Display raincloudplots - two figures in poster
[cb] = cbrewer('qual','Set3',12,'pchip');
cl(1,:) = cb(4,:);
cl(2,:) = cb(1,:);
fig_position = [0 0 1 1]; % coordinates for figures

if I_roi_distr ~= 0
    % Figure 5 A and B in poster
    f7 = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
    ax1 = subplot(1,3,1)
    h1= raincloud_plot('X', tsnr_img{1}(I_mask) , 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
    h2= raincloud_plot('X', tsnr_img{2}(I_mask), 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
    h3= raincloud_plot('X', tsnr_img{3}(I_mask), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
    h4= raincloud_plot('X', tsnr_img{4}(I_mask), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0);
    legend([h1{1} h2{1} h3{1} h4{1}], {'Echo 2', 'pre-T2{\ast}', 'pre-tSNR', 'pre-TE'})
    title('tSNR distributions within whole brain mask')
    set(gca,'XLim', [-10 350], 'YLim', [-.02 .02]);
    hold(ax1, 'on')
    plot(ax1, [0 0],[0 -0.2],'k')
    hold(ax1, 'off')
    box off

    ax2 = subplot(1,3,2)
    h1= raincloud_plot('X', tsnr_img{1}(I_roi_distr) , 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
    h2= raincloud_plot('X', tsnr_img{2}(I_roi_distr), 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
    h3= raincloud_plot('X', tsnr_img{3}(I_roi_distr), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
    h4= raincloud_plot('X', tsnr_img{4}(I_roi_distr), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0);
    legend([h1{1} h2{1} h3{1} h4{1}], {'Echo 2', 'pre-T2{\ast}', 'pre-tSNR', 'pre-TE'})
    title(['tSNR distributions within ' roi_text] )
    set(gca,'XLim', [-10 300], 'YLim', [-.02 .02]);
    hold(ax2, 'on')
    plot(ax2, [0 0],[0 -0.2],'k')
    hold(ax2, 'off')
    box off

    ax3 = subplot(1,3,3)
    h1= raincloud_plot('X', perc_diff{1}(I_roi_distr) , 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
    h2= raincloud_plot('X', perc_diff{2}(I_roi_distr), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
    h3= raincloud_plot('X', perc_diff{3}(I_roi_distr), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
    legend([h1{1} h2{1} h3{1}], {'pre-T2{\ast}', 'pre-tSNR', 'pre-TE'})
    title(['Percentage difference distributions within ' roi_text] )
    set(gca,'XLim', [-50 percentage_threshold], 'YLim', [-.04 .04]);
    hold(ax3, 'on')
    plot(ax3, [0 0],[0 -0.4],'k')
    hold(ax3, 'off')
    box off
else
% Figure 5 A and B in poster
    f7 = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
    ax1 = subplot(1,1,1)
    h1= raincloud_plot('X', tsnr_img{1}(I_mask) , 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
    h2= raincloud_plot('X', tsnr_img{2}(I_mask), 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
    h3= raincloud_plot('X', tsnr_img{3}(I_mask), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
    h4= raincloud_plot('X', tsnr_img{4}(I_mask), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0);
    legend([h1{1} h2{1} h3{1} h4{1}], {'Echo 2', 'pre-T2{\ast}', 'pre-tSNR', 'pre-TE'})
    title('tSNR distributions within whole brain mask')
    set(gca,'XLim', [-10 350], 'YLim', [-.02 .02]);
    hold(ax1, 'on')
    plot(ax1, [0 0],[0 -0.2],'k')
    hold(ax1, 'off')
    box off

end


if ~exist(distr_saveAs_fn, 'file')
    print(f7, distr_saveAs_fn,'-dpng', '-r0');
end
