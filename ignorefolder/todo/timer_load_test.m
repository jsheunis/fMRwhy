f = figure;
global load;
load = struct;
load.str1 = '';
load.str2 = '.';
load.str3 = '..';
load.str4 = '...';

global txt_preproc1;
txt_preproc1 = uicontrol('Parent', f,...
    'Style','text',...
    'String','-',...
    'Units', 'Normalized',...
    'Position',[0.5 0.5 0.2 0.3],...
    'HorizontalAlignment', 'left',...
    'fontsize', 30);

% uicontrol('Style','text', 'String',char(hex2dec('2713')), ...
%     'Units','normalized', 'Position',[0 0 1 1], ...
%     'FontSize',30)
% 
% uicontrol('Style','text', 'String',char(hex2dec('2713')), ...
%     'Units','normalized', 'Position',[0 0 1 1], ...
%     'FontSize',30)

period = 0.5;

t = timer;
t.UserData = 1;
t.TimerFcn = @updateStrings;
t.Period = 0.5;
t.ExecutionMode = 'fixedRate';

start(t)

function updateStrings(t,~)
global load txt_preproc1;

if t.UserData == 5
    t.UserData = 1;
end

txt_preproc1.String = load.(['str' num2str(t.UserData)]);


t.UserData = t.UserData + 1;
drawnow;

end