function fmrwhy_util_thePlotROI(functional_fn, mask_fn, roi_fn, task_info, trace_info, saveAs_fn)
%
% This is a Matlab script that uses custom code and spm12 routines to plot
% a version of THE PLOT from Jonathan Power. See:
% - https://www.jonathanpower.net/2017-ni-the-plot.html
% - https://www.ncbi.nlm.nih.gov/pubmed/27510328
% (code has been partly adopted from Power's script)
%
% THE PLOT is a useful visual quality checking tool for fMRI timeseries
% data. It is a 2D plot of voxel intensities over time, with voxels ordered
% in bins (GM, WM, CSF) to indicate some directionality in the data. Motion
% and physiological noise are more easily visually inspected with this
% tool, especially if viewed next to other quality traces, like FD, DVARS
% or physiological recordings.

% 
% UPDATES:
% - 23 July 2019: added voxel ordering options according to RO and GSO
% methods in: https://www.biorxiv.org/content/10.1101/662726v1.full

% ------------------------------------------------------------------------ %
% functional_fn = 4D timeseries; typically smoothed (NOT realigned nor slice time corrected)
% mask_fn = binary brain mask; exclude voxels outside of this mask
% roi_fn = filename of ROI
% task_info = info about task paradigm for plotting purposes
% trace_info = traces for for plotting purposes
% ------------------------------------------------------------------------ %


% ThePlot settings:
intensity_scale = [-6 6]; % scaling for plot image intensity, see what works
fontsizeL = 15;
fontsizeM = 13;
visibility = 'off';

%ordering = 0; % ordering of voxels in ThePlotSPM;
% 0 = Random order via standard Matlab indexing (RO)
% 1 = Grey matter signal ordering (GSO) (according to the Pearsons correlation between voxel timeseries and global signal )
% 2 = Cluster-similarity ordering (CO) - NOT IMPLEMENTED YET

% roi = ''; % figure generated for voxels in supplied ROI - NOT IMPLEMENTED YET

% ------ %
% STEP 1 %
% ------ %
% Load timeseries data, get image and time dimensions, transform data
nii = nii_tool('load', functional_fn);
data_4D = double(nii.img);
[Ni, Nj, Nk, Nt] = size(data_4D); % [Ni x Nj x Nk x Nt]
data_2D = reshape(data_4D, Ni*Nj*Nk, Nt); %[voxels, time]
data_2D = data_2D'; %[time, voxels]

% Load mask data
mask_nii = nii_tool('load', mask_fn);
mask_3D = mask_nii.img; % [Ni x Nj x Nk]
mask_2D = reshape(mask_3D, Ni*Nj*Nk, 1); % [Ni*Nj*Nk x 1]
mask_I = find(mask_2D); % [Nmaskvoxels x 1]
mask_I = mask_I'; % [1 x Nmaskvoxels]

% Load ROI data
roi_nii = nii_tool('load', roi_fn);
roi_3D = roi_nii.img; % [Ni x Nj x Nk]
roi_3D = fmrwhy_util_createBinaryImg(roi_3D, 0.1);
roi_2D = reshape(roi_3D, Ni*Nj*Nk, 1); % [Ni*Nj*Nk x 1]
roi_I = find(roi_2D); % [Nroivoxels x 1]
roi_I = roi_I'; % [1 x Nroivoxels]

% TODO: make sure that roi_I and mask_I are the correct format, length, orientation, etc.

% Load task data
[task_time_course, convolved_ttc] = fmrwhy_util_createBlockParadigm(Nt, task_info.TR, task_info.onsets, task_info.durations, task_info.precision);

% ------ %
% STEP 2 %
% ------ %

% Remove linear and quadratic trend per voxel
data_2D_detrended = fmrwhy_util_detrend(data_2D, 2); %[time, voxels]

% Calculate mean
data_2D_mean = nanmean(data_2D_detrended); %[1, voxels]

% Demean
data_2D_demeaned = data_2D_detrended - data_2D_mean; %[time, voxels]

% Calculate standard deviation
data_2D_stddev = std(data_2D_detrended); %[1, voxels]

% Calculate percentage signal change: [I(t) - mean(I)]/mean(I)*100
data_2D_psc = 100*(data_2D_detrended./data_2D_mean) - 100;
data_2D_psc(isnan(data_2D_psc)) = 0;
F_2D_psc = data_2D_psc';


% ------ %
% STEP 3 %
% ------ %

% Prepare timeseries
% TODO: note that global signal gets assigned the colour code of grey matter in plot, these are not the same. update this.
scale_f = 0.8;
%gs = fmrwhy_util_detrend(confounds_struct.global_signal, 1);
%gs = 6+fmrwhy_util_scale(gs,-scale_f,scale_f);
%wm = fmrwhy_util_detrend(confounds_struct.white_matter, 1);
%wm = 3.5+fmrwhy_util_scale(wm,-scale_f,scale_f);
%csf = fmrwhy_util_detrend(confounds_struct.csf, 1);
%csf = 1+fmrwhy_util_scale(csf,-scale_f,scale_f);
%%wm = 2+zscore(confounds_struct.white_matter);
%%wm = fmrwhy_util_detrend(wm, 1);
%%csf = 1+zscore(confounds_struct.csf);
%%csf = fmrwhy_util_detrend(csf, 1);
%
%physmat_fn = fullfile(options.func_dir_qc, ['PhysIO_task-' task '_run-' run], 'physio.mat');
%physio = load(physmat_fn);
%card = physio.physio.ons_secs.c;
%resp = 2+physio.physio.ons_secs.r;
%%resp = fmrwhy_util_scale(resp,-0.5,0.5);
%t_phys = 1:numel(card);

% ------ %
% STEP 4 %
% ------ %

% Get screen size for plotting
scr_size = get(0,'ScreenSize');
dist = scr_size(4);
if scr_size(3) < dist
    dist = scr_size(3);
end

% Specify colors
colors_hex = {'#F25F5C', '#FFE066', '#247BA0', '#2E294E', '#70C1B3', '#6B2737'};
for i = 1:numel(colors_hex)
    colors_rgb{i} = sscanf(colors_hex{i}(2:end),'%2x%2x%2x',[1 3])/255;
end

% Order voxels
ROI_signals = F_2D_psc(roi_I, :);
ROI_time_series = mean(ROI_signals, 1);
Rcorr = corr(ROI_signals', ROI_time_series');
[R_sorted, I_sorted] = sort(Rcorr,'descend');
all_img = ROI_signals(I_sorted, :);
% Plot
f = figure('units','normalized','outerposition',[0 0 1 1], 'visible', visibility);
% Ax1 - The Plot
ax1 = subplot(7,1,[4:7]);
imagesc(ax1, all_img); colormap(gray); caxis(intensity_scale);
xlabel(ax1, 'fMRI volumes','fontsize',fontsizeM);
xlim(ax1,[0 Nt]);
set(ax1,'Yticklabel',[]);
fmrwhy_util_stretchAx(ax1);
axpos = ax1.Position;
% Ax2 - Task paradigm
ax2 = subplot(7,1,1);
im_task = imagesc(ax2, task_time_course'); colormap(gray);
%plot(ax2, gs, 'LineWidth', 2, 'Color', colors_rgb{1});
%hold(ax2, 'on')
%plot(ax2, wm, 'LineWidth', 2, 'Color', colors_rgb{2});
%plot(ax2, csf, 'LineWidth', 2, 'Color', colors_rgb{3});
%hold(ax2, 'off')
xlim(ax2,[0 Nt]);
%ylim(ax2,[0 7.5])
set(ax2,'Xticklabel',[]);
set(ax2,'Yticklabel',[]);
ax2.XAxis.Visible = 'off';
ax2.YAxis.Visible = 'off';
fmrwhy_util_stretchAx(ax2);
ax2.Position(1) = axpos(1); ax2.Position(3) = axpos(3);

% Ax3 - Convolved task
ax3 = subplot(7,1,2);
plot(ax3, convolved_ttc, 'LineWidth', 2, 'Color', colors_rgb{1});
%hold(ax3, 'on')
%plot(ax3, card, 'LineWidth', 2, 'Color', colors_rgb{5});
%hold(ax3, 'off')
xlim(ax3,[0 Nt]);
ylim(ax3,[-0.5 1.5]);
set(ax3,'Xticklabel',[]);
set(ax3,'Yticklabel',[]);
ax3.XAxis.Visible = 'off';
ax3.YAxis.Visible = 'off';
fmrwhy_util_stretchAx(ax3)
ax3.Position(1) = axpos(1); ax3.Position(3) = axpos(3);
% Ax4 - ROI signal
ax4 = subplot(7,1,3);
ROI_ts = fmrwhy_util_scale(ROI_time_series,-scale_f,scale_f);
plot(ax4, ROI_time_series, 'LineWidth', 2, 'Color', colors_rgb{3});
xlim(ax4,[0 Nt]);
ylim(ax4,[-2 2]);
set(ax4,'Xticklabel',[]);
set(ax4,'Yticklabel',[]);
ax4.XAxis.Visible = 'off';
ax4.YAxis.Visible = 'off';
fmrwhy_util_stretchAx(ax4);
ax4.Position(1) = axpos(1); ax4.Position(3) = axpos(3);
%% [left bottom width height]
ax1pos = ax1.OuterPosition;
ax2pos = ax2.OuterPosition;
ax3pos = ax3.OuterPosition;
ax4pos = ax4.OuterPosition;
ax3.Position(2) = ax2pos(2) - ax3pos(4);
ax4.Position(2) = ax2pos(2) - ax3pos(4) - ax4pos(4);
ax1.Position(4) = ax2pos(2) - ax3pos(4) - ax4pos(4) - ax1pos(2);
%txt_legend = {'Global signal', 'WM signal', 'CSF signal', 'Respiration', 'Heart rate', 'Framewise displacement'};
txt_legend = {'Global signal', 'WM signal', 'CSF signal', 'Respiration', 'Heart rate', 'Framewise displacement'};
xt = [2 2 2 2000 2000 2]; yt = [7.2 4.7 2.2 3.4 1.3 1];

h_txt(1) = text(ax3, 2, 1.2, 'Task', 'Color', colors_rgb{1},'FontSize',14);
h_txt(2) = text(ax4, 2, 1.5, 'ROI signal', 'Color', colors_rgb{3},'FontSize',14);

% Save figure
if ~exist(saveAs_fn, 'file')
    print(f, saveAs_fn,'-dpng')
end