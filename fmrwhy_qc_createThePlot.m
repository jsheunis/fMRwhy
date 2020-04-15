function fmrwhy_qc_createThePlot(bids_dir, sub, ses, task, run, echo, options)
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

% Load/create required defaults
% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% Update workflow params with subject functional derivative filenames
options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

% ThePlot settings:
intensity_scale = options.theplot.intensity_scale; % scaling for plot image intensity, see what works
fontsizeL = 15;
fontsizeM = 13;

%ordering = 0; % ordering of voxels in ThePlotSPM;
% 0 = Random order via standard Matlab indexing (RO)
% 1 = Grey matter signal ordering (GSO) (according to the Pearsons correlation between voxel timeseries and global signal )
% 2 = Cluster-similarity ordering (CO) - NOT IMPLEMENTED YET
% roi = ''; % figure generated for voxels in supplied ROI - NOT IMPLEMENTED YET

% Timeseries to use for ThePlot
functional4D_fn = options.sfunctional_fn;

% Get image information from functional data
nii = nii_tool('load', functional4D_fn);
data_4D = nii.img;
[Ni, Nj, Nk, Nt] = size(data_4D);

% Load multiple confound regressors
confounds_struct = tdfread(options.confounds_fn);
confounds_mat = struct2array(confounds_struct);

% Load mask data
masks = fmrwhy_util_loadMasks(bids_dir, sub);

% Calculate stats data (including BOLD percentage signal change)
stats = fmrwhy_qc_calculateStats(bids_dir, sub, functional4D_fn);

% Load data
F_2D_psc = reshape(stats.data_4D_psc, Ni*Nj*Nk, Nt); %[voxels, time]
I_GM = masks.GM_mask_I;
I_WM = masks.WM_mask_I;
I_CSF = masks.CSF_mask_I;
I_brain = masks.brain_mask_I;

% Prepare timeseries
% TODO: note that global signal gets assigned the colour code of grey matter in plot, these are not the same. update this.
scale_f = 0.8;
gs = fmrwhy_util_detrend(confounds_struct.global_signal, 1);
gs = 6+fmrwhy_util_scale(gs,-scale_f,scale_f);
wm = fmrwhy_util_detrend(confounds_struct.white_matter, 1);
wm = 3.5+fmrwhy_util_scale(wm,-scale_f,scale_f);
csf = fmrwhy_util_detrend(confounds_struct.csf, 1);
csf = 1+fmrwhy_util_scale(csf,-scale_f,scale_f);
%wm = 2+zscore(confounds_struct.white_matter);
%wm = fmrwhy_util_detrend(wm, 1);
%csf = 1+zscore(confounds_struct.csf);
%csf = fmrwhy_util_detrend(csf, 1);

physmat_fn = fullfile(options.func_dir_qc, ['PhysIO_task-' task '_run-' run], 'physio.mat');
physio = load(physmat_fn);
card = physio.physio.ons_secs.c;
resp = 2+physio.physio.ons_secs.r;
%resp = fmrwhy_util_scale(resp,-0.5,0.5);
t_phys = 1:numel(card);

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

% ---
% Figure 0 - ordering 0
% ---
% Order voxels
order_text = '(RO)';
GM_img = F_2D_psc(I_GM, :);
WM_img = F_2D_psc(I_WM, :);
CSF_img = F_2D_psc(I_CSF, :);
all_img = [GM_img; WM_img; CSF_img];
% Plot
f = figure('units','normalized','outerposition',[0 0 1 1]);
% Ax1 - The Plot
ax1 = subplot(7,1,[4:7]);
imagesc(ax1, all_img); colormap(gray); caxis(intensity_scale);
xlabel(ax1, 'fMRI volumes','fontsize',fontsizeM)
% Ax1 - Add patches if ordered by tissue compartment
hold on;
gm_pos = numel(I_GM);
gwm_pos = numel(I_GM) + numel(I_WM);
gwmcsf_pos = numel(I_GM) + numel(I_WM) + numel(I_CSF);
x_gm = [0 0 1 1]; y_gm = [0 gm_pos gm_pos 0];
x_wm = [0 0 1 1]; y_wm = [gm_pos gwm_pos gwm_pos gm_pos];
x_csf = [0 0 1 1]; y_csf = [gwm_pos gwmcsf_pos gwmcsf_pos gwm_pos];
p1 = patch(x_gm, y_gm, colors_rgb{1},'EdgeColor','none');
p2 = patch(x_wm, y_wm, colors_rgb{2},'EdgeColor','none');
p3 = patch(x_csf, y_csf, colors_rgb{3},'EdgeColor','none');
hold off;
xlim(ax1,[0 Nt])
set(ax1,'Yticklabel',[]);
fmrwhy_util_stretchAx(ax1)
axpos = ax1.Position;
% Ax2 - Brain signals, etc
ax2 = subplot(7,1,1);
plot(ax2, gs, 'LineWidth', 2, 'Color', colors_rgb{1});
hold(ax2, 'on')
plot(ax2, wm, 'LineWidth', 2, 'Color', colors_rgb{2});
plot(ax2, csf, 'LineWidth', 2, 'Color', colors_rgb{3});
hold(ax2, 'off')
xlim(ax2,[0 Nt])
ylim(ax2,[0 7.5])
set(ax2,'Xticklabel',[]);
set(ax2,'Yticklabel',[]);
ax2.XAxis.Visible = 'off';
ax2.YAxis.Visible = 'off';
fmrwhy_util_stretchAx(ax2)
ax2.Position(1) = axpos(1); ax2.Position(3) = axpos(3);
% Ax3 - Physiology
ax3 = subplot(7,1,2);
plot(ax3, resp, 'LineWidth', 2, 'Color', colors_rgb{4});
hold(ax3, 'on')
plot(ax3, card, 'LineWidth', 2, 'Color', colors_rgb{5});
hold(ax3, 'off')
xlim(ax3,[0 numel(card)])
ylim(ax3,[-1 4])
set(ax3,'Xticklabel',[]);
set(ax3,'Yticklabel',[]);
ax3.XAxis.Visible = 'off';
ax3.YAxis.Visible = 'off';
fmrwhy_util_stretchAx(ax3)
ax3.Position(1) = axpos(1); ax3.Position(3) = axpos(3);
% Ax4 - Framewise displacement
ax4 = subplot(7,1,3);
plot(ax4, confounds_struct.framewise_displacement, 'LineWidth', 2, 'Color', colors_rgb{6});
xlim(ax4,[0 Nt])
ylim(ax4,[0 1])
set(ax4,'Xticklabel',[]);
set(ax4,'Yticklabel',[]);
ax4.XAxis.Visible = 'off';
ax4.YAxis.Visible = 'off';
fmrwhy_util_stretchAx(ax4)
ax4.Position(1) = axpos(1); ax4.Position(3) = axpos(3);
%% [left bottom width height]
ax1pos = ax1.OuterPosition;
ax2pos = ax2.OuterPosition;
ax3pos = ax3.OuterPosition;
ax4pos = ax4.OuterPosition;
ax3.Position(2) = ax2pos(2) - ax3pos(4);
ax4.Position(2) = ax2pos(2) - ax3pos(4) - ax4pos(4);
ax1.Position(4) = ax2pos(2) - ax3pos(4) - ax4pos(4) - ax1pos(2);
txt_legend = {'Global signal', 'WM signal', 'CSF signal', 'Respiration', 'Heart rate', 'Framewise displacement'};
xt = [2 2 2 2000 2000 2]; yt = [7.2 4.7 2.2 3.4 1.3 0.5];
for i = 1:6
    if i<4
        h_txt(i) = text(ax2, xt(i), yt(i), txt_legend{i}, 'Color', colors_rgb{i},'FontSize',14);
    elseif i>3 && i<6
        h_txt(i) = text(ax3, xt(i), yt(i), txt_legend{i}, 'Color', colors_rgb{i},'FontSize',14);
    else
        h_txt(i) = text(ax4, xt(i), yt(i), txt_legend{i}, 'Color', colors_rgb{i},'FontSize',14);
    end
end
% Save figure
theplot_fn = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-RO_grayplot.png']);
%if ~exist(theplot_fn, 'file')
    print(f,theplot_fn,'-dpng')
%end

% ---
% Figure 1 - ordering 1
% ---
% Order voxels
GM_signals = F_2D_psc(I_GM, :);
GM_time_series = mean(GM_signals, 1);
F_masked_psc = F_2D_psc(I_brain, :);
Rcorr = corr(F_masked_psc', GM_time_series');
[R_sorted, I_sorted] = sort(Rcorr,'descend');
all_img = F_masked_psc(I_sorted, :);
% Plot
f = figure('units','normalized','outerposition',[0 0 1 1]);
% Ax1 - The Plot
ax1 = subplot(7,1,[4:7]);
imagesc(ax1, all_img); colormap(gray); caxis(intensity_scale);
xlabel(ax1, 'fMRI volumes','fontsize',fontsizeM)
xlim(ax1,[0 Nt])
set(ax1,'Yticklabel',[]);
fmrwhy_util_stretchAx(ax1)
axpos = ax1.Position;
% Ax2 - Brain signals, etc
ax2 = subplot(7,1,1);
plot(ax2, gs, 'LineWidth', 2, 'Color', colors_rgb{1});
hold(ax2, 'on')
plot(ax2, wm, 'LineWidth', 2, 'Color', colors_rgb{2});
plot(ax2, csf, 'LineWidth', 2, 'Color', colors_rgb{3});
hold(ax2, 'off')
xlim(ax2,[0 Nt])
ylim(ax2,[0 7.5])
set(ax2,'Xticklabel',[]);
set(ax2,'Yticklabel',[]);
ax2.XAxis.Visible = 'off';
ax2.YAxis.Visible = 'off';
fmrwhy_util_stretchAx(ax2)
ax2.Position(1) = axpos(1); ax2.Position(3) = axpos(3);
% Ax3 - Physiology
ax3 = subplot(7,1,2);
plot(ax3, resp, 'LineWidth', 2, 'Color', colors_rgb{4});
hold(ax3, 'on')
plot(ax3, card, 'LineWidth', 2, 'Color', colors_rgb{5});
hold(ax3, 'off')
xlim(ax3,[0 numel(card)])
ylim(ax3,[-1 4])
set(ax3,'Xticklabel',[]);
set(ax3,'Yticklabel',[]);
ax3.XAxis.Visible = 'off';
ax3.YAxis.Visible = 'off';
fmrwhy_util_stretchAx(ax3)
ax3.Position(1) = axpos(1); ax3.Position(3) = axpos(3);
% Ax4 - Framewise displacement
ax4 = subplot(7,1,3);
plot(ax4, confounds_struct.framewise_displacement, 'LineWidth', 2, 'Color', colors_rgb{6});
xlim(ax4,[0 Nt])
ylim(ax4,[0 1])
set(ax4,'Xticklabel',[]);
set(ax4,'Yticklabel',[]);
ax4.XAxis.Visible = 'off';
ax4.YAxis.Visible = 'off';
fmrwhy_util_stretchAx(ax4)
ax4.Position(1) = axpos(1); ax4.Position(3) = axpos(3);
%% [left bottom width height]
ax1pos = ax1.OuterPosition;
ax2pos = ax2.OuterPosition;
ax3pos = ax3.OuterPosition;
ax4pos = ax4.OuterPosition;
ax3.Position(2) = ax2pos(2) - ax3pos(4);
ax4.Position(2) = ax2pos(2) - ax3pos(4) - ax4pos(4);
ax1.Position(4) = ax2pos(2) - ax3pos(4) - ax4pos(4) - ax1pos(2);
txt_legend = {'Global signal', 'WM signal', 'CSF signal', 'Respiration', 'Heart rate', 'Framewise displacement'};
xt = [2 2 2 2000 2000 2]; yt = [7.2 4.7 2.2 3.4 1.3 0.5];
for i = 1:6
    if i<4
        h_txt(i) = text(ax2, xt(i), yt(i), txt_legend{i}, 'Color', colors_rgb{i},'FontSize',14);
    elseif i>3 && i<6
        h_txt(i) = text(ax3, xt(i), yt(i), txt_legend{i}, 'Color', colors_rgb{i},'FontSize',14);
    else
        h_txt(i) = text(ax4, xt(i), yt(i), txt_legend{i}, 'Color', colors_rgb{i},'FontSize',14);
    end
end
% Save figure
theplot_fn = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-GSO_grayplot.png']);
%if ~exist(theplot_fn, 'file')
    print(f,theplot_fn,'-dpng')
%end



% Order 2:
        % TODO: implement this correctly
%        % From: https://gist.github.com/benfulcher/f243bdc4ab80a7351f083f35bdd4db63
%        dataMatrix = F_masked_psc';
%        distanceMetric = 'Euclidean';
%        linkageMethod = 'average';
%%         dataMatrixNorm = zscore(dataMatrix); % normalize columns
%        dataMatrixNorm = dataMatrix;
%        R = pdist(dataMatrixNorm,distanceMetric); % Pairwise distances
%        links = linkage(R,linkageMethod); % Do hierarchical linkage
%        f = figure('color','w');
%        set(gcf,'Visible','off'); % suppress figure output
%        [~,~,ord] = dendrogram(links,0);
%        close(f); % close the invisible figure used for the dendrogram
%        dataMatrixClustered = dataMatrixNorm(ord,:); % Reorder rows
%        all_img = dataMatrixClustered';