
%--------------------------------------------------------------------------
% SETUP OF VISUALISATIONS
%--------------------------------------------------------------------------

% Create figure for real-time visualisation
fig_plot = figure('units','normalized','outerposition',[0 0 1 0.66]);
Fplot = {};
for e = 1:Ne
    Fplot{e} = reshape(F{e}, Nx, Ny, Nz, Ndyn);
    rFplot{e} = reshape(rF{e}, Nx, Ny, Nz, Ndyn);

end
slice_plot = 26;

ax1 = subplot(1,3,1);
ax2 = subplot(1,3,2);
ax3 = subplot(1,3,3);

im1 = imagesc(ax1, rot90(squeeze(Fplot{1}(:,:,slice_plot, 84)),rotateVal)); colormap gray;
im2 = imagesc(ax2, rot90(squeeze(rFplot{1}(:,:,slice_plot, 84)),rotateVal)); colormap gray;
im3 = imagesc(ax3, rot90(squeeze(Fplot{3}(:,:,slice_plot, 1)),rotateVal)); colormap gray;
t1 = text(ax1, 6, 6, num2str(1),'FontSize',16,'Color','white');
drawnow;

%for i=1:Ndyn
%    i
%    new1 = rot90(squeeze(Fplot{1}(:,:,slice_plot, i)),rotateVal);
%    new2 = rot90(squeeze(Fplot{2}(:,:,slice_plot, i)),rotateVal);
%    new3 = rot90(squeeze(Fplot{3}(:,:,slice_plot, i)),rotateVal);
%    set(im1, 'CData', new1);
%    set(im2, 'CData', new2);
%    set(im3, 'CData', new3);
%    set(t1,'String',num2str(i))
%    drawnow;
%    pause(0.1)
%end