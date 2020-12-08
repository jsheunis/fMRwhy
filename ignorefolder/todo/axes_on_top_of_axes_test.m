

fig = figure;

ax1 = axes('Parent', fig,...
    'Visible', 'off',...
    'Position',[.2 .28 .6 .62],...
    'Color', 'k',...
    'FontSize', 18);

ax2 = axes('Parent', fig,...
    'Visible', 'off',...
    'Position',[.2 .28 .6 .62],...
    'Color', 'k',...
    'FontSize', 18);
%%
x = 1:10;
y = x.^2;
x2 = 1:10;
y2 = -x.^1.5;
p1 = plot(ax1, x, y);
%%
p2 = plot(ax2, x2, y2);

%%
set(p2, 'Xdata', [1 2.5 3.3 4.5 5.4 6 7 8 9 10]);
drawnow;
% p2 = plot(ax2, x2, y2);

%%

ax1.Visible = 'off';
p1.Visible = 'off';

%%
ax2.Visible = 'off';
p2.Visible = 'off';