function fmrwhy_util_plotTimeseriesMovie(func1, func2, func3, slice_plot)
    % Create figure for real-time visualisation
    % fig_plot = figure;
    fig_plot = figure('units', 'normalized', 'outerposition', [0 0 1 0.66]);
    rotateVal = 1;

    Ndyn = numel(func1);
    ax1 = subplot(1, 3, 1);
    ax2 = subplot(1, 3, 2);
    ax3 = subplot(1, 3, 3);

    im1 = imagesc(ax1, rot90(squeeze(func1{1}(slice_plot, :, :)), rotateVal));
    colormap gray;
    im2 = imagesc(ax2, rot90(squeeze(func2{1}(slice_plot, :, :)), rotateVal));
    colormap gray;
    im3 = imagesc(ax3, rot90(squeeze(func3{1}(slice_plot, :, :)), rotateVal));
    colormap gray;
    t1 = text(ax1, 6, 6, num2str(1), 'FontSize', 16, 'Color', 'white');
    drawnow;

    for i = 1:Ndyn
        new1 = rot90(squeeze(func1{i}(slice_plot, :, :)), rotateVal);
        %    new2 = rot90(squeeze(func2{i}(slice_plot,:,:)),rotateVal);
        %    new3 = rot90(squeeze(func3{i}(slice_plot,:,:)),rotateVal);
        set(im1, 'CData', new1);
        %    set(im2, 'CData', new2);
        %    set(im3, 'CData', new3);
        set(t1, 'String', num2str(i));
        drawnow;
        pause(0.02);
    end
