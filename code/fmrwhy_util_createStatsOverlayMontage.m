function output = fmrwhy_util_createStatsOverlayMontage(template_img, stats_img, roi_img, columns, rotate, str, clrmp, visibility, shape, cxs, stats_clrmp, roi_rgbcolors, clrbar, saveAs_fn)
% Function to create montages of images overlaid on a template image

% Structure to save output
output = struct;

% Prepare plotting options
alpha = 0.1;
plot_contour = 1;
if roi_rgbcolors == 0
    roi_rgbcolors = [255,255,191; 215,25,28; 253,174,97; 171,217,233; 44,123,182]/255;
    %roi_rgbcolors = [215,25,28; 253,174,97; 255,255,191; 171,217,233; 44,123,182]/255;
else
    roi_rgbcolors = roi_rgbcolors/255;
end
% Get screen size for plotting
scr_size = get(0,'ScreenSize');
dist = scr_size(4);
if scr_size(3) < dist
    dist = scr_size(3);
end

% Create background montage
montage_template = fmrwhy_util_createMontage(template_img, columns, rotate, 'Template volume', clrmp, 'off', shape, cxs);
% Create stats montage
if ~isempty(stats_img)
    montage_stats = fmrwhy_util_createMontage(stats_img, columns, rotate, 'Overlay', stats_clrmp, 'off', shape, 'auto');
end
% Create ROI montage(s)
if ~isempty(roi_img)
    montage_rois = {};
    if iscell(roi_img)
        for i=1:numel(roi_img)
            montage_rois{i} = fmrwhy_util_createMontage(roi_img{i}, columns, rotate, 'Overlay', clrmp, 'off', shape, 'auto');
        end
    else
        montage_rois{1} = fmrwhy_util_createMontage(roi_img, columns, rotate, 'Overlay', clrmp, 'off', shape, 'auto');
    end
end
% Prepare colours for ROIs
[Nimx, Nimy] = size(montage_template.whole_img);
oo = ones(Nimx, Nimy);
zz = zeros(Nimx, Nimy);
%red = cat(3, oo, zz, zz);
%green = cat(3, zz, oo, zz);
%blue = cat(3, zz, oo, oo);

% Create figure with background montage and overlaid motages (stats and rois)
f = figure('units','normalized','outerposition',[0 0 1 1], 'visible', visibility);
ax1 = axes('Parent', f);
im1 = imagesc(ax1, montage_template.whole_img);
colormap(ax1, clrmp);
if ~isempty(cxs)
    caxis(ax1, cxs);
end
ax1 = fmrwhy_util_stretchAx(ax1);
ax1 = fmrwhy_util_removeTicksAx(ax1);

ax2 = axes('Parent',f);

if ~isempty(stats_img)
    stats = montage_stats.whole_img;
    imagesc(ax2, stats, 'alphadata', stats>0);
    colormap(ax2, stats_clrmp);
else
    stats = nan(size(montage_template.whole_img));
    imagesc(ax2, stats,'alphadata', stats>0);
    colormap(ax2, clrmp);
end

set(ax2, 'color','none','visible','off');
ax2 = fmrwhy_util_stretchAx(ax2);
linkaxes([ax1 ax2]);


% Plot ROIs
if ~isempty(roi_img)
    hold(ax2, 'on');
    for i=1:numel(montage_rois)
        rbgclr = roi_rgbcolors(i,:);
        clr = cat(3, rbgclr(1)*oo, rbgclr(2)*oo, rbgclr(3)*oo);
        imC = imagesc(ax2, clr);
        set(imC, 'AlphaData', alpha*montage_rois{i}.whole_img);
        if plot_contour
            bound_whole_bin = bwboundaries(montage_rois{i}.whole_img);
            Nblobs_bin = numel(bound_whole_bin);
            for b = 1:Nblobs_bin
                p = plot(ax2, bound_whole_bin{b,1}(:,2), bound_whole_bin{b,1}(:,1), 'color', rbgclr, 'LineWidth', 1.5);
            end
        end
    end
    hold(ax2, 'off');
end

% Add custom colorbar for stat map
if clrbar
    [s1, s2] = size(montage_template.whole_img);
    x1 = s2 - s2/columns/2 - round(s2/columns/8);
    x2 = s2 - s2/columns/2 + round(s2/columns/8);
    x = [x1 x2 x2 x1];
    y1 = s1 - s2/columns/4;
    y2 = y1 - s2/columns/2;
    y = [y1 y1 y2 y2];

    if ~isempty(cxs)
        cmin = cxs(1);
        cmax = cxs(2);
    else
        cmin = min(min(stats));
        cmax = max(max(stats));
    end
    C = [cmin cmin cmax cmax];
    pp = patch(ax2, x,y,C);
    pp.LineWidth = 1;
    pp.EdgeColor = [91, 92, 91]/255;
    % Add t-values to colorbar
    t1 = text(ax1, x2+2, y1, num2str(round(cmin, 2)),'FontSize',13,'Color','white');
    t2 = text(ax1, x2+2, y2, num2str(round(cmax, 2)),'FontSize',13,'Color','white');
end


% Save and close figures
output.f = f;
if saveAs_fn ~= 0
    print(f, saveAs_fn,'-dpng', '-r0')
end
% Close necessary figure handles
close(montage_template.f)
if ~isempty(stats_img)
    close(montage_stats.f)
end
if ~isempty(roi_img)
    for i=1:numel(montage_rois)
        close(montage_rois{i}.f)
    end
end
if strcmp(visibility, 'off')
    close(f);
end
