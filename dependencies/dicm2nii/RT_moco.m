function RT_moco()
% Display and save motion information at real time. It also shows images and the
% progress of EPI/DTI while scanning, and allows to check motion information for
% previous runs/patients.
%
% To make this work, you will need:
%  1. Set up shared folder for ../incoming_DICOM at the computer running RT_moco.
%  2. Set up real time image transfer at Siemens console.
% 
% Here are some example command at Siemens console to set up real time transfer:
%  net use Y: \\yourIP\yourPath\incoming_DICOM yourPassWd /USER:yourUserName
%  xedit -f C:\MedCom\config\Ice\IceConfig.evp -n ICE.CONFIG.OnlineSendConfiguration -p OnlineTargetPort -v -1
%  xedit -f C:\MedCom\config\Ice\IceConfig.evp -n ICE.CONFIG.OnlineSendConfiguration -p OnlineTargetPath -v Y:\
%  xedit -f C:\MedCom\config\Ice\IceConfig.evp -n ICE.CONFIG.OnlineSendConfiguration -p OnlineSendIMA -v 1
% 
% Limitations of current version:
%  1. Tested only for Siemens Prisma at CCBBI at OSU. 
%  2. Requires later matlab (2017a is the earliest version tested).
%  3. Ideal for portrait display monitor setup for nows.

% 200207 xiangrui.li at gmail.com first working version inspired by FIRMM

if ~exist('./log/', 'dir'), mkdir('./log/'); end % folder to save subj.mat

% Create/re-use GUI and start timer.
fh = findall(0, 'Type', 'figure', 'Tag', 'RT_moco');
if ~isempty(fh)
    figure(fh); hs = guidata(fh);
    if hs.timer.Running == "off", start(hs.timer); end % in case it stopped
    return;
end
res = get(0, 'ScreenSize');
fh = figure('mc'*[256 1]'); clf(fh);
set(fh, 'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off', ... 
	'DockControls', 'off', 'CloseRequestFcn', @closeFig, 'Color', 'w', ...
    'Name', 'Real Time Image Monitor.', 'Position', [2 60 1076 res(4)-120], ...
    'Tag', 'RT_moco', 'UserData', struct('FD', {{}}, 'DV', {{}}, 'hdr', {{}}));
hs.fig = fh;

h = uimenu(fh, 'Label', '&Patient');
hs.menu(1) = uimenu(h, 'Label', 'Load Patient', 'Callback', @loadSubj);
hs.menu(2) = uimenu(h, 'Label', 'Redo Patient', 'Callback', @redoSubj);
hs.menu(3) = uimenu(h, 'Label', 'Close Patient', 'Callback', @closeSubj);

h = uimenu(fh, 'Label', '&Series');
uimenu(h, 'Label', 'View Selected Series in 3D', 'Callback', @view_3D);
uimenu(h, 'Label', 'Overlay Selected Series onto Anatomy', 'Callback', @overlay);
hs.derived = uimenu(h, 'Label', 'Skip DERIVED Series', 'Checked', 'on', ...
    'Callback', @toggleChecked, 'Separator', 'on');
hs.SBRef = uimenu(h, 'Label', 'Skip *_SBRef Series', 'Callback', @toggleChecked, 'Checked', 'on');

h = uimenu(fh, 'Label', '&View');
uimenu(h, 'Label', 'Increase Brightness', 'Callback', @setCLim);
uimenu(h, 'Label', 'Decrease Brightness', 'Callback', @setCLim);
uimenu(h, 'Label', 'Show FD plot', 'Callback', @toggleFD, 'Separator', 'on');
h = uimenu(h, 'Label', '&FD Threshold');
for i = [0.1:0.1:0.5 0.8 1 2], uimenu(h, 'Label', num2str(i), 'Callback', @FD_yLim); end

dy = 0.12 * (0:3);
hs.ax = axes(fh, 'Units', 'normalized', 'Position', [0.05 0.46 0.9 0.18], ...
    'NextPlot', 'add', 'XLim', [0.5 300.5], 'UserData', dy, ...
    'TickDir', 'out', 'TickLength', 0.002*[1 1], 'ColorOrder', [0 0 1; 1 0 1]);
xlabel(hs.ax, 'Instance Number');
hs.slider = uicontrol(fh, 'Units', 'normalized', 'Position', [0.035 0.645 0.93 0.018], ...
    'Style', 'slider', 'Value', 1, 'Min', 1, 'Max', 300, 'Callback', @sliderCB, ...
    'BackgroundColor', 0.5*[1 1 1], 'SliderStep', [1 1]./300);

yyaxis left; ylabel(hs.ax, 'DVARS');
set(hs.ax, 'YTick', dy, 'YLim', dy([1 4]));
c3 = [0 0.8 0;  0.8 0.8 0;  0.5 0 0];
for i = 3:-1:1
    rectangle(hs.ax, 'Position', [0.5 dy(i) 2000 dy(i+1)-dy(i)], ...
        'FaceColor', c3(i,:), 'EdgeColor', c3(i,:), 'LineWidth', 0.01);
end
hs.dv = plot(hs.ax, 0, '.:');

yyaxis right; ylabel(hs.ax, 'Framewise Displacement (mm)');
set(hs.ax, 'YTick', 0:0.4:1.2, 'YLim', [0 1.2]);

txt = @(a)text(hs.ax, 'Units', 'normalized', 'Position', a, 'FontSize', 12, ...
    'BackgroundColor', 'w', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
hs.pct(1) = txt([0.995 0.32]); hs.pct(2) = txt([0.995 0.65]);

hs.fd = plot(hs.ax, 0, '.:', 'Visible', 'off');
hs.ax.YAxis(2).Visible = 'off';
vars = {'Description' 'Series' 'Instances' '<font color="#00cc00">Green</font>' ...
    '<font color="#cccc00">Yellow</font>' 'MeanFD'};
hs.table = uitable(fh, 'Units', 'normalized', 'Position', [0.05 0.01 0.9 0.4], ...
    'FontSize', 14, 'RowName', [], 'ColumnWidth', {370 100 130 100 100 120}, ...
    'ColumnName', strcat('<html><h2>', vars, '</h2></html>'), ...
    'CellSelectionCallback', @tableCB);

ax = axes(fh, 'Units', 'normalized', 'Position', [0.7 0.67 0.29 0.33], 'Visible', 'off');
hs.subj = text(ax, 'Position', [0.85 1], 'FontSize', 24, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Interpreter', 'none');
hs.series = text(ax, 'Position', [0.85 0], 'FontSize', 18, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Interpreter', 'none');

ax = axes(fh, 'Units', 'normalized', 'Position', [0.05 0.67 0.65 0.33], ...
    'YDir', 'reverse', 'Visible', 'off');
hs.img = image(ax, 'CData', inf(2), 'CDataMapping', 'scaled');
axis equal; colormap gray;
hs.instnc = text(ax, 'Units', 'normalized', 'Position', [0.99 0.01], 'Color', 'y', ...
    'FontSize', 14, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

hs.timer = timer('StartDelay', 0.2, 'ObjectVisibility', 'off', ...
    'StopFcn', @saveResult, 'TimerFcn', @doSeries, 'UserData', fh);
fh.HandleVisibility = 'callback'; % fh.Resize = 'off'; 
guidata(fh, hs);
start(hs.timer);

%% TimerFunc: do a series if avail, then call stopFunc to save result
function doSeries(obj, ~)
hs = guidata(obj.UserData);
hs.fig.Name(end) = 78 - hs.fig.Name(end); % indication of timer call
[f, iRun] = next_series(hs);
if isempty(f), return; end
set(hs.menu, 'Enable', 'off'); set([hs.table hs.slider], 'Enable', 'inactive');

asc_header = dicm2nii('', 'asc_header', 'func_handle');
dict = dicm_dict('', {'Rows' 'Columns' 'BitsAllocated' 'InstanceNumber'});
rob = java.awt.Robot(); key = java.awt.event.KeyEvent.VK_SHIFT;
rob.keyPress(key); rob.keyRelease(key); % wake up screen

hs.series.UserData = iRun; % needed for next_series()
s = dicm_hdr_wait(sprintf('%s/001_%06.f_000001.dcm', f, iRun));
if isempty(s), return; end % non-image dicom, skip series
if iRun == 1 % first series: reset GUI
    closeSubj(hs.fig);
    hs.subj.String = strrep(s.PatientName, ' ', '_'); hs.subj.UserData = f;
    close(findall(0, 'Type', 'figure', 'Tag', 'nii_viewer')); % last subj if any
end

if hs.derived.Checked=="on" && contains(s.ImageType, 'DERIVED'), return; end
if hs.SBRef.Checked=="on" && endsWith(s.SeriesDescription, '_SBRef'), return; end

hs.series.String = seriesInfo(s);
try 
    nSL = s.CSAImageHeaderInfo.NumberOfImagesInMosaic; % EPI | DTI mosaic
catch % T1, T2, fieldmap etc: show info/img only
    nSL = asc_header(s, 'sSliceArray.lSize'); % 2D
    if nSL==1, nSL = asc_header(s, 'sKSpace.lImagesPerSlab'); end % 3D
    iSL = ceil(nSL/2); % try middle slice for better view
    nam = sprintf('%s/001_%06.f_%06.f.dcm', f, iRun, iSL);
    nIN = numel(dir([nam(1:end-9) '*.dcm']));
    if nIN<nSL, pause(nSL/50); nIN = numel(dir([nam(1:end-9) '*.dcm'])); end
    if now-getfield(dir(s.Filename), 'datenum') < 1/1440 % within 1 min
        nIN = nIN * asc_header(s, 'sSliceArray.lConc'); % 2D T2
    end
    if ~exist(nam, 'file')
        iSL = 1; nam = sprintf('%s/001_%06.f_%06.f.dcm', f, iRun, iSL);
    end
    init_series(hs, s, nIN);
    set_img(hs.img, dicm_img(nam));
    hs.slider.Value = iSL;
    hs.instnc.String = num2str(iSL);
    return;
end

isDTI = contains(s.ImageType, '\DIFFUSION');
if isDTI, nIN = asc_header(s, 'sDiffusion.lDiffDirections') + 1; % free-dir too
else, nIN = asc_header(s, 'lRepetitions') + 1;
end
if isempty(nIN) || endsWith(s.SeriesDescription, '_SBRef'), nIN = 1; end
mos = dicm_img(s); s.PixelData = mos;
img = mos2vol(mos, nSL);
p = refVol(img, [s.PixelSpacing' s.SpacingBetweenSlices]);

ijk = round(p.R0 \ p.mm + 1); % informative voxels
ind = sub2ind(size(img), ijk(1,:), ijk(2,:), ijk(3,:));
img0 = double(mos);
mn = mean(img0(ind));
s.CLim = mn + std(img0(ind)*3); if isDTI, s.CLim = s.CLim / 2; end
set_img(hs.img, mos, s.CLim);
init_series(hs, s, nIN);
viewer = findall(0, 'Type', 'figure', 'Tag', 'nii_viewer');
if nIN<6 && isempty(viewer), overlay(hs.fig); end

R1 = inv(p.R0);
m6 = zeros(2,6);
hs.fd.YData(2:end) = nan; hs.dv.YData(2:end) = nan;

nextSeries = sprintf('%s/001_%06.f_000001.dcm', f, iRun+1);
try dt = s.RepetitionTime/1000 + 6; catch, dt = 9; end % for stopped series
for i = 2:nIN
    nam = sprintf('%s/001_%06.f_%06.f.dcm', f, iRun, i);
    tEnd = now + dt/86400;
    while ~exist(nam, 'file')
        if ~isempty(dir(nextSeries)) || now>tEnd, return; end
        pause(0.2); 
    end
    s = dicm_hdr_wait(nam, dict); iN = s.InstanceNumber;
    mos = dicm_img(s);
    img = mos2vol(mos, nSL);
    hs.img.CData = mos; hs.instnc.String = num2str(iN);
    hs.slider.Value = iN; % show progress
    if isDTI, hs.table.Data{1,3} = i; hs.dv.YData(iN) = 0; continue; end
    p.F.Values = smooth_mc(img, p.sz);
    [m6(2,:), R1] = moco_estim(p, R1);
    a = abs(m6(2,:) - m6(1,:)); m6(1,:) = m6(2,:);
    hs.fd.YData(iN) = sum([a(1:3) a(4:6)*50]); % 50mm: head radius
    
    img = double(mos);
    a = img(ind) - img0(ind); % use only edge voxles: faster and more sensitive
    hs.dv.YData(iN) = sqrt(a*a' / numel(a)) / mn;
    img0 = img;
    a = hs.dv.YData(1:iN); a = a(~isnan(a));
    dy = hs.ax.UserData; fd = hs.fd.YData(1:iN);
    N = {numel(a) sum(a<dy(2)) sum(a<dy(3))};
    hs.table.Data(1,3:6) = [N mean(fd(~isnan(fd)))];
    for j = 1:2, hs.pct(j).String = sprintf('%.3g%%', N{j+1}/N{1}*100); end
    if iN>=nIN, return; end % ISSS alike
    % drawnow; % update instance for offline test
end

%% Reshape mosaic into volume, remove padded zeros
function vol = mos2vol(mos, nSL)
nMos = ceil(sqrt(nSL)); % nMos x nMos tiles for Siemens
[nr, nc] = size(mos); % number of row & col in mosaic
nr = nr / nMos; nc = nc / nMos; % number of row and col in slice
vol = zeros([nr nc nSL], class(mos));
for i = 1:nSL
    % r =    mod(i-1, nMos) * nr + (1:nr); % 2nd slice is tile(2,1)
    % c = floor((i-1)/nMos) * nc + (1:nc);
    r = floor((i-1)/nMos) * nr + (1:nr); % 2nd slice is tile(1,2)
    c =    mod(i-1, nMos) * nc + (1:nc);
    vol(:, :, i) = mos(r, c);
end

%% Initialize GUI for a new series
function init_series(hs, s, nIN)
hs.dv.YData = zeros(nIN,1); hs.fd.YData = zeros(nIN,1);
hs.table.Data = [{s.SeriesDescription s.SeriesNumber nIN 0 0 0}; hs.table.Data];
set(hs.slider, 'Max', nIN, 'Value', 1, 'UserData', s.Filename(1:end-10));
if nIN==1, hs.slider.Visible = 'off';
else, set(hs.slider, 'SliderStep', [1 1]./(nIN-1)); hs.slider.Visible = 'on';
end
hs.ax.XLim(2) = nIN + 0.5; 
hs.fig.UserData.hdr{end+1} = s; % 1st instance with CLim and maybe image
set([hs.instnc hs.pct], 'String', '');
figure(hs.fig); drawnow; % bring GUI front if needed

%% Set img and img axis
function set_img(hImg, img, CLim)
if nargin<3 || isempty(CLim) || isnan(CLim)
    img0 = double(img(img>100));
    if isempty(img0), img0 = 1; end
    CLim = mean(img0) + std(img0)*2;
end
d = size(img) + 0.5;
set(hImg.Parent, 'CLim', [0 CLim], 'XLim', [0.5 d(2)], 'YLim', [0.5 d(1)]);
hImg.CData = img;

%% get some series information
function c = seriesInfo(s)
c{1} = s.SeriesDescription;
c{2} = sprintf('Series %g', s.SeriesNumber);
c{3} = datestr(datenum(s.AcquisitionTime, 'HHMMSS.fff'), 'HH:MM:SS AM');
try c{4} = sprintf('TR = %g', s.RepetitionTime); catch, end

%% toggle FD display on/off
function toggleFD(h, ~)
hs = guidata(h);
if hs.fd.Visible == "on"
    set([hs.fd hs.ax.YAxis(2)], 'Visible', 'off');
    h.Label = 'Show FD plot';
else
    set([hs.fd hs.ax.YAxis(2)], 'Visible', 'on');
    h.Label = 'Hide FD plot';
end

%% Set FD plat y-axis limit
function FD_yLim(h, ~)
hs = guidata(h);
dy = str2double(h.Label) * (0:3);
yyaxis(hs.ax, 'right'); set(hs.ax, 'YTick', dy, 'YLim', dy([1 4]));

%% Table-click callback: show moco/series info and image if avail
function tableCB(h, evt)
if isempty(evt.Indices) || evt.Indices(1,2)>2, return; end
hs = guidata(h);
C = h.Data;
i = evt.Indices(1,1);
iR = size(C,1) - i + 1;
hs.fd.YData = hs.fig.UserData.FD{iR};
hs.dv.YData = hs.fig.UserData.DV{iR};
for j = 1:2, hs.pct(j).String = sprintf('%.3g%%', C{i,j+3}/C{i,3}*100); end

hs.instnc.String = '';
hs.series.String = C{i,1}; % in case hdr not saved
try s = hs.fig.UserData.hdr{iR}; catch, set_img(hs.img, inf(2), 1); return; end

nIN = sum(~isnan(hs.dv.YData));
iIN = ceil(nIN/2); % start with middle Instance if avail
nam = sprintf('%s%06g.dcm', s.Filename(1:end-10), iIN);
if ~exist(nam, 'file'), iIN = 1; nam = s; end
hs.instnc.String = num2str(iIN);
set(hs.slider, 'Max', nIN, 'Value', iIN, 'UserData', s.Filename(1:end-10));
if nIN == 1, hs.slider.Visible = 'off';
else, set(hs.slider, 'SliderStep', [1 1]./(nIN-1)); hs.slider.Visible = 'on';
end
hs.ax.XLim(2) = numel(hs.dv.YData) + 0.5;
hs.series.String = seriesInfo(s);
try CLim = s.CLim; catch, CLim = []; end
set_img(hs.img, dicm_img(nam), CLim);

%% Load subj data to review
function loadSubj(h, ~)
[fname, pName] = uigetfile('./log/*.mat', 'Select a Subject to load');
if isnumeric(fname), return; end
load([pName '/' fname], 'T3');
hs = guidata(h);
hs.fig.UserData = T3.Properties.UserData; 
DV = hs.fig.UserData.DV;
N = size(T3, 1);
C = flip(table2cell(T3), 1); C(:,6) = C(:,3);
dy = hs.ax.UserData;
for i = 1:N
    a = DV{N-i+1}; a = a(~isnan(a));
    C(i,3:5) = {numel(a) sum(a<dy(2)) sum(a<dy(3))};
end 
hs.table.Data = C;
s = hs.fig.UserData.hdr{end};
hs.subj.String = strrep(s.PatientName, ' ', '_');
hs.subj.UserData = fileparts(s.Filename);
hs.series.UserData = str2double(s.Filename(end+(-16:-11)));
tableCB(hs.table, struct('Indices', [1 1])); % show top series

%% close subj
function closeSubj(h, ~)
hs = guidata(h);
hs.table.Data = {};
hs.img.CData = inf(2);
hs.fd.YData = 0; hs.dv.YData = 0;
hs.subj.UserData = '';
hs.fig.UserData = struct('FD', {{}}, 'DV', {{}}, 'hdr', {{}});
set([hs.subj hs.series hs.instnc hs.pct], 'String', '');

%% Re-do current subj: useful in case of error during a session
function redoSubj(h, ~)
hs = guidata(h);
subj = hs.subj.String;
if isempty(subj), return; end
if ~exist(hs.subj.UserData, 'dir')
    fprintf(2, 'Image for %s deleted?\n', subj);
    return;
end
try delete(['./log/' subj '*.mat']); catch, end
hs.table.Data = {}; % quick visual sign

%% Get reference vol info. Adapted from nii_moco.m
function p = refVol(img, pixdim)
d = size(img);
p.R0 = diag([pixdim 1]); % no need for real xform_mat here
p.R0(1:3, 4) = -pixdim .* (d/2); % make center voxel [0 0 0]

sz = pixdim;
if all(abs(diff(sz)/sz(1))<0.05) && sz(1)>2 && sz(1)<4 % 6~12mm
    sz = 3; % iso-voxel, 2~4mm res, simple fast smooth
else
    sz = 9 ./ sz'; % 9 mm seems good
end

% resample ref vol to isovoxel (often lower-res)
d0 = d-1;
dd = 4 ./ pixdim; % use 4 mm grid for alignmen
[i, j, k] = ndgrid(0:dd(1):d0(1)-0.5, 0:dd(2):d0(2)-0.5, 0:dd(3):d0(3)-0.5);
I = [i(:) j(:) k(:)]';
a = rng('default'); I = I + rand(size(I))*0.5; rng(a); % used by spm
V = smooth_mc(img, sz);
F = griddedInterpolant({0:d0(1), 0:d0(2), 0:d0(3)}, V, 'linear', 'none');
V0 = F(I(1,:), I(2,:), I(3,:)); % ref: 1 by nVox
I(4,:) = 1; % 0-based ijk: 4 by nVox
I = p.R0 * I; % xyz of ref voxels

% compute derivative to each motion parameter in ref vol
dG = zeros(6, numel(V0));
dd = 1e-6; % delta of motion parameter, value won't affect dG much
R0i = inv(p.R0); % speed up a little
for i = 1:6
    p6 = zeros(6,1); p6(i) = dd; % change only 1 of 6
    J = R0i * rigid_mat(p6) * I; %#ok<*MINV>
    dG(i,:) = F(J(1,:), J(2,:), J(3,:)) - V0; % diff now
end
dG = dG / dd; % derivative

% choose voxels with larger derivative for alignment: much faster
a = sum(dG.^2); % 6 derivatives has similar range
ind = a > std(a(~isnan(a)))/10; % arbituray threshold. Also exclude NaN
p.dG = dG(:, ind);
p.V0 = V0(ind);
p.mm = I(:, ind);
F.GridVectors = {0:d(1)-1, 0:d(2)-1, 0:d(3)-1};
p.F = F;
p.sz = sz;

%% motion correction to ref-vol. From nii_moco.m
function [m6, rst] = moco_estim(p, R)
mss0 = inf;
rst = R;
for iter = 1:64
    J = R * p.mm; % R_rst*J -> R0*ijk
    V = p.F(J(1,:), J(2,:), J(3,:));
    ind = ~isnan(V); % NaN means out of range
    dV = p.V0(ind) - V(ind);
    mss = dV*dV' / numel(dV); % mean(dV.^2)
    if mss > mss0, break; end % give up and use previous R
    rst = R; % accecpt only if improving
    if 1-mss/mss0 < 1e-6, break; end % little effect, stop
    
    a = p.dG(:, ind);
    p6 = (a * a') \ (a * dV'); % dG(:,ind)'\dV' estimate p6 from current R
    R = R * rigid_mat(p6); % inv(inv(rigid_mat(p6)) * inv(R_rst))
    mss0 = mss;
end

R = p.R0 * rst; % inv(R_rst / Rref)
m6 = -[R(1:3, 4)' atan2(R(2,3), R(3,3)) asin(R(1,3)) atan2(R(1,2), R(1,1))];

%% Translation (mm) and rotation (deg) to 4x4 R. Order: ZYXT
function R = rigid_mat(p6)
ca = cosd(p6(4:6)); sa = sind(p6(4:6));
rx = [1 0 0; 0 ca(1) -sa(1); 0 sa(1) ca(1)]; % 3D rotation
ry = [ca(2) 0 sa(2); 0 1 0; -sa(2) 0 ca(2)];
rz = [ca(3) -sa(3) 0; sa(3) ca(3) 0; 0 0 1];
R = rx * ry * rz;
R = [R p6(1:3); 0 0 0 1];

%% Simple gaussian smooth for motion correction, sz in unit of voxels
function out = smooth_mc(in, sz)
out = double(in);
if all(abs(diff(sz)/sz(1))<0.05) && abs(sz(1)-round(sz(1)))<0.05 ...
        && mod(round(sz(1)),2)==1
    out = smooth3(out, 'gaussian', round(sz)); % sz odd integer
    return; % save time for special case
end

d = size(in);
I = {1:d(1) 1:d(2) 1:d(3)};
n = sz/3;
if numel(n)==1, n = n*[1 1 1]; end
J = {1:n(1):d(1) 1:n(2):d(2) 1:n(3):d(3)};
intp = 'linear';
F = griddedInterpolant(I, out, intp);
out = smooth3(F(J), 'gaussian'); % sz=3
F = griddedInterpolant(J, out, intp);
out = F(I);

%% new series or new subj: result saved as ./log/subj.mat
% The subject folders (yyyymmdd.SubjName.SubjID) is under ../incoming_DICOM/
% The dcm file names from Siemens push are always in format of
% 001_000001_000001.dcm. The first num is always 1 (StudyID?), and both second
% (series) and third (instance) are always continuous.
function [f, iRun] = next_series(hs)
if ~isempty(hs.subj.UserData) % check new run for current subj
    f = hs.subj.UserData;
    iRun = hs.series.UserData + 1;
    if ~isempty(dir(sprintf('%s/001_%06.f_000001.dcm', f, iRun))), return; end
end
rootDir = '../incoming_DICOM/';
dirs = dir([rootDir '20*']); % check new subj
dirs(~[dirs.isdir]) = [];
f = ''; iRun = 1;
for i = numel(dirs):-1:1
    subj = regexp(dirs(i).name, '(?<=20\d{6}\.).*?(?=\.)', 'match', 'once');
    if exist(['./log/' subj '.mat'], 'file'), continue; end
    f = [rootDir dirs(i).name]; return;
end

if ispc || mod(now,1) > hs.timer.StartDelay*2/86400; return; end % mid-night?
for i = 1:numel(dirs)
    if now-dirs(i).datenum < 2, continue; end % keep data for 2 days
    try rmdir([rootDir dirs(i).name], 's');
    catch me, disp(me.message);
    end
end

%% Wait till the file is copied completely
function s = dicm_hdr_wait(varargin)
tEnd = now + 1/86400; % wait up to 1 second
while 1
    s = dicm_hdr(varargin{:});
    try if s.PixelData.Start+s.PixelData.Bytes <= s.FileSize, return; end; end %#ok
    if now>tEnd, s = []; return; end % give up
    pause(0.1);
end

%% User closing GUI: stop and delete timer
function closeFig(fh, ~)
hs = guidata(fh);
try hs.timer.StopFcn = ''; stop(hs.timer); delete(hs.timer); catch, end
delete(fh);

%% menu callback for both DERIVED and _SBRef
function toggleChecked(h, ~)
if h.Checked == "on", h.Checked = "off"; else, h.Checked = "on"; end

%% Increase/Decrease image CLim
function setCLim(h, ~)
hs = guidata(h);
ax = hs.img.Parent;
if startsWith(h.Label, 'Increase'), ax.CLim(2) = ax.CLim(2)*0.8;
else, ax.CLim(2) = ax.CLim(2)*1.2;
end

%% show series in nii_viewer
function view_3D(h, ~)
hs = guidata(h);
nams = dir([hs.slider.UserData '*.dcm']);
if isempty(nams), return; end
nams = strcat(nams(1).folder, '/', {nams.name});
nii = dicm2nii(nams, ' ', 'no_save');
nii_viewer(nii);

%% overlay series onto T1w or Scout if avail
function overlay(h, ~)
hs = guidata(h);
hdrs = hs.fig.UserData.hdr;
if isempty(hdrs), return; end
is3D = @(c)c.MRAcquisitionType=="3D" && ~contains(c.ImageType, 'DERIVED');
is3D = cellfun(is3D, hdrs);
if ~any(is3D), is3D = cellfun(@(c)c.MRAcquisitionType=="3D", hdrs); end
if sum(is3D)>1
    isT1 = cellfun(@(c)contains(c.SequenceName, 'fl3d1'), hdrs);
    isT1 = isT1 & is3D;
    if any(isT1), is3D = isT1; end
end
if ~any(is3D), view_3D(h); return; end % no T1, just show in nii_viewer
is3D = find(is3D, 1, 'last');
nams = dir([hdrs{is3D}.Filename(1:end-10) '*.dcm']);
nams = strcat(nams(1).folder, '/', {nams.name});
T1w = dicm2nii(nams, ' ', 'no_save');
nams = dir([hs.slider.UserData '*.dcm']);
nams = strcat(nams(1).folder, '/', {nams.name});
epi = dicm2nii(nams, ' ', 'no_save');
nii_viewer(T1w, epi);
fh = findall(0, 'type', 'figure', 'tag', 'nii_viewer');
nii_viewer('LocalFunc', 'nii_viewer_cb', [], [], 'center', fh(1));

%% slider callback: show img if avail
function sliderCB(h, ~)
hs = guidata(h);
if isempty(h.UserData), return; end
dict = dicm_dict('', {'SamplesPerPixel' 'Rows' 'Columns' 'BitsAllocated' 'InstanceNumber'});
h.Value = round(h.Value);
s = dicm_hdr(sprintf('%s%06.f.dcm', h.UserData, h.Value), dict);
if isempty(s), return; end
hs.img.CData = dicm_img(s);
hs.instnc.String = num2str(s.InstanceNumber);

%% Timer StopFunc: Save result, start timer with delay, even after error.
function saveResult(obj, ~)
hs = guidata(obj.UserData);
set([hs.menu hs.table hs.slider], 'Enable', 'on');
if size(hs.table.Data,1) > numel(hs.fig.UserData.FD) % new series to save?
    hs.fig.UserData.FD{end+1} = hs.fd.YData;
    hs.fig.UserData.DV{end+1} = hs.dv.YData;
    T3 = cell2table(flip(hs.table.Data(:,[1 2 6]), 1), ...
        'VariableNames', {'Description' 'SeriesNumber' 'MeanFD'});
    T3.Properties.UserData = hs.fig.UserData;
    save(['./log/' hs.subj.String], 'T3');
end

if isempty(next_series(hs)), obj.StartDelay = 5; else, obj.StartDelay = 0.2; end
start(obj);

%%