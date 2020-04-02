function output = fmrwhy_util_createOverlayMontage(template_img, overlay_img, columns, rotate, str, clrmp, visibility, shape, cxs, saveAs_fn)
% Function to create montages of images/rois overlaid on a template image

% Structure to save output
output = struct;
alpha = 0.1;
plot_contour = 1;
rgbcolors = [255,255,191; 215,25,28; 253,174,97; 171,217,233; 44,123,182]/255;
%rgbcolors = [215,25,28; 253,174,97; 255,255,191; 171,217,233; 44,123,182]/255;

% Create background montage
montage_template = fmrwhy_util_createMontage(template_img, columns, rotate, 'Template volume', clrmp, 'off', shape, cxs);
%caxis(montage_template.ax, cxs);
% Get screen size for plotting
scr_size = get(0,'ScreenSize');
dist = scr_size(4);
if scr_size(3) < dist
    dist = scr_size(3);
end

% Create figures with background montage and overlaid masks
%f(i) = figure('units','pixels','outerposition',[0 0 dist dist]);
f = figure('units','normalized','outerposition',[0 0 1 1], 'visible', visibility);
im1 = imagesc(montage_template.whole_img);
colormap(clrmp);
colorbar;
caxis(cxs);
ax = gca;
%if cxs ~= 0
%    caxis(ax, cxs);
%end
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
hold(ax, 'on')
[Nimx, Nimy] = size(montage_template.whole_img);
oo = ones(Nimx, Nimy);
zz = zeros(Nimx, Nimy);
red = cat(3, oo, zz, zz);
green = cat(3, zz, oo, zz);
blue = cat(3, zz, oo, oo);

if iscell(overlay_img)
    for i=1:numel(overlay_img)
        montage_overlay{i} = fmrwhy_util_createMontage(overlay_img{i}, columns, rotate, 'Overlay', clrmp, 'off', shape, 'auto');
    end
else
    montage_overlay = {};
    montage_overlay{1} = fmrwhy_util_createMontage(overlay_img, columns, rotate, 'Overlay', clrmp, 'off', shape, 'auto');
end


for i=1:numel(montage_overlay)
    rbgclr = rgbcolors(i,:);
    clr = cat(3, rbgclr(1)*oo, rbgclr(2)*oo, rbgclr(3)*oo);
    imC = imagesc(ax, clr);
    set(imC, 'AlphaData', alpha*montage_overlay{i}.whole_img);
    if plot_contour
        bound_whole_bin = bwboundaries(montage_overlay{i}.whole_img);
        Nblobs_bin = numel(bound_whole_bin);
        for b = 1:Nblobs_bin
        p = plot(ax, bound_whole_bin{b,1}(:,2), bound_whole_bin{b,1}(:,1), 'color', rbgclr, 'LineWidth', 1);
        end
    end
end



hold(ax, 'off');
set(ax,'xtick',[])
set(ax,'xticklabel',[])
set(ax,'ytick',[])
set(ax,'yticklabel',[])
set(ax,'ztick',[])
set(ax,'zticklabel',[])

output.ax = ax;
output.f = f;

if saveAs_fn ~= 0
    print(f, saveAs_fn,'-dpng', '-r0')
end
% Close necessary figure handles
close(montage_template.f)
for i=1:numel(montage_overlay)
    close(montage_overlay{i}.f)
end
if strcmp(visibility, 'off')
    close(f);
end
