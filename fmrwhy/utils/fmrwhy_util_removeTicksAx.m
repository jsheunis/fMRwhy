function ax = fmrwhy_util_removeTicksAx(ax)
    % Remove all axes ticks
    set(ax, 'xtick', []);
    set(ax, 'xticklabel', []);
    set(ax, 'ytick', []);
    set(ax, 'yticklabel', []);
    set(ax, 'ztick', []);
    set(ax, 'zticklabel', []);
