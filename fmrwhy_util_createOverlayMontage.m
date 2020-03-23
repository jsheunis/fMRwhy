function output = fmrwhy_util_createOverlayMontage(template_img, overlay_img, columns, rotate, str, clrmp, visibility, shape)
% Function to create montages of images/rois overlaid on a template image

% Structure to save output
output = struct;

% Create background montage
montage_EPI = fmrwhy_util_createMontage(template_img, columns, rotate, 'Template volume', clrmp, 'off', shape);

% Get screen size for plotting
scr_size = get(0,'ScreenSize');
dist = scr_size(4);
if scr_size(3) < dist
    dist = scr_size(3);
end

% Create figures with background montage and overlaid masks
%f(i) = figure('units','pixels','outerposition',[0 0 dist dist]);
f = figure('units','normalized','outerposition',[0 0 1 1]);
im1 = imagesc(montage_EPI.whole_img);
colormap('gray');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
hold(ax, 'on')
[Nimx, Nimy] = size(montage_EPI.whole_img);
oo = ones(Nimx, Nimy);
zz = zeros(Nimx, Nimy);
red = cat(3, oo, zz, zz);
montage_overlay = fmrwhy_util_createMontage(overlay_img, columns, rotate, 'Overlay', clrmp, 'off', shape);
bound_whole_bin = bwboundaries(montage_overlay.whole_img);
Nblobs_bin = numel(bound_whole_bin);
for b = 1:Nblobs_bin
p = plot(ax, bound_whole_bin{b,1}(:,2), bound_whole_bin{b,1}(:,1), 'r', 'LineWidth', 1);
end
imC = imagesc(ax, red);
set(imC, 'AlphaData', 0.2*montage_overlay.whole_img);
hold(ax, 'off');
set(ax,'xtick',[])
set(ax,'xticklabel',[])
set(ax,'ytick',[])
set(ax,'yticklabel',[])
set(ax,'ztick',[])
set(ax,'zticklabel',[])
print(f, 'testy_overlay','-dpng', '-r0')

