function fmrwhy_util_computeROImeasures(functional_fn, roi_img, task_info, tsnr_fn, title_str, saveAs_fn, options)

% Functional_fn = filename for slice time corrected, realigned and smoothed data.

% Compute task region info:
%    - Signal dropout/intensity comparison (histogram?)
%    - The plot (without CSF / WM), with convolved task paradigm
%    - Task region trace (spatial average over time) with convolved task paradigm
%    - tSNR average comparison and dotplot / raincloud
%    - Percentage signal change comparison (from voxel timeseries or average)

%/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-rrightMotor_roi.nii

% ThePlot settings:
intensity_scale = options.theplot.intensity_scale; % scaling for plot image intensity, see what works
fontsizeL = 15;
fontsizeM = 13;
%
%% Timeseries to use for ThePlot
%functional4D_fn = options.sfunctional_fn;

% Get image information
func_spm = spm_vol(functional_fn);
Nt = numel(func_spm);
Ni = func_spm(1).dim(1);
Nj = func_spm(1).dim(2);
Nk = func_spm(1).dim(3);

% Load multiple confound regressors
%confounds_struct = tdfread(options.confounds_fn);
%confounds_mat = struct2array(confounds_struct);

% Load task data
[task_time_course, convolved_ttc] = fmrwhy_util_createBlockParadigm(Nt, task_info.TR, task_info.onsets, task_info.durations, task_info.precision);

% Load mask data
roi_2D = reshape(roi_img, 1, Ni*Nj*Nk); %[1, voxels]
I_roi = find(roi_2D);
Nvox_roi = numel(I_roi)

% Load tSNR data
% (this tsnr image provided as argument will typically be calculated from slice time corrected and realigned data)
tsnr = spm_read_vols(spm_vol(tsnr_fn));

% Load functional timeseries
data_4D = spm_read_vols(func_spm);
data_2D = reshape(data_4D, Ni*Nj*Nk, Nt); %[voxels, time]
data_2D = data_2D'; %[time, voxels]

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

% Prepare timeseries
% TODO: note that global signal gets assigned the colour code of grey matter in plot, these are not the same. update this.
scale_f = 0.8;
%gs = fmrwhy_util_detrend(confounds_struct.global_signal, 1);
%gs = 6+fmrwhy_util_scale(gs,-scale_f,scale_f);
%wm = fmrwhy_util_detrend(confounds_struct.white_matter, 1);
%wm = 3.5+fmrwhy_util_scale(wm,-scale_f,scale_f);
%csf = fmrwhy_util_detrend(confounds_struct.csf, 1);
%csf = 1+fmrwhy_util_scale(csf,-scale_f,scale_f);



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
% Figure 1 - ordering 1
% ---
% Order voxels
ROI_signals = F_2D_psc(I_roi, :);
ROI_time_series = mean(ROI_signals, 1);
Rcorr = corr(ROI_signals', ROI_time_series');
[R_sorted, I_sorted] = sort(Rcorr,'descend');
all_img = ROI_signals(I_sorted, :);
% Plot
f = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'on');
% Ax1 - The Plot
ax1 = subplot(7,1,[4:7]);
imagesc(ax1, all_img); colormap(gray); caxis(intensity_scale);
xlabel(ax1, 'fMRI volumes','fontsize',fontsizeM)
xlim(ax1,[0 Nt])
set(ax1,'Yticklabel',[]);
fmrwhy_util_stretchAx(ax1)
axpos = ax1.Position;
% Ax2 - Task paradigm
ax2 = subplot(7,1,1);
im_task = imagesc(ax2, task_time_course'); colormap(gray);
%plot(ax2, gs, 'LineWidth', 2, 'Color', colors_rgb{1});
%hold(ax2, 'on')
%plot(ax2, wm, 'LineWidth', 2, 'Color', colors_rgb{2});
%plot(ax2, csf, 'LineWidth', 2, 'Color', colors_rgb{3});
%hold(ax2, 'off')
xlim(ax2,[0 Nt])
%ylim(ax2,[0 7.5])
set(ax2,'Xticklabel',[]);
set(ax2,'Yticklabel',[]);
ax2.XAxis.Visible = 'off';
ax2.YAxis.Visible = 'off';
fmrwhy_util_stretchAx(ax2)
ax2.Position(1) = axpos(1); ax2.Position(3) = axpos(3);
title(title_str)
% Ax3 - Convolved task
ax3 = subplot(7,1,2);
plot(ax3, convolved_ttc, 'LineWidth', 2, 'Color', colors_rgb{1});
%hold(ax3, 'on')
%plot(ax3, card, 'LineWidth', 2, 'Color', colors_rgb{5});
%hold(ax3, 'off')
xlim(ax3,[0 Nt])
ylim(ax3,[-0.5 1.5])
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
xlim(ax4,[0 Nt])
ylim(ax4,[-2 2])
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
%txt_legend = {'Global signal', 'WM signal', 'CSF signal', 'Respiration', 'Heart rate', 'Framewise displacement'};
txt_legend = {'Global signal', 'WM signal', 'CSF signal', 'Respiration', 'Heart rate', 'Framewise displacement'};
xt = [2 2 2 2000 2000 2]; yt = [7.2 4.7 2.2 3.4 1.3 1];

h_txt(1) = text(ax3, 2, 1.2, 'Task', 'Color', colors_rgb{1},'FontSize',14);
h_txt(2) = text(ax4, 2, 1.5, 'ROI signal', 'Color', colors_rgb{3},'FontSize',14);
%for i = 1:6
%    if i<4
%        h_txt(i) = text(ax2, xt(i), yt(i), txt_legend{i}, 'Color', colors_rgb{i},'FontSize',14);
%    elseif i>3 && i<6
%        h_txt(i) = text(ax3, xt(i), yt(i), txt_legend{i}, 'Color', colors_rgb{i},'FontSize',14);
%    else
%
%    end
%end
% Save figure
%theplot_fn = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-GSO_grayplot.png']);
if saveAs_fn ~= 0
    if ~exist(saveAs_fn, 'file')
        print(f,saveAs_fn,'-dpng')
    end
end



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