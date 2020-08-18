
roi = 1; % left motor
signal = 1; % echo 2

figure;
hold on;
for i = 1:6
    whole_brain_ts =  cell2mat(signals_raw(i,:));
    roi_voxels_ts = whole_brain_ts(I_roi{roi},:);
    roi_ts = mean(roi_voxels_ts, 1);
    roi_ts_z = zscore(roi_ts);
    plot(roi_ts_z);
end

legend(gca,{'Echo 2','Pre-tSNR','Pre-T2star','Pre-TE','RT-T2star','T2star FIT'}, 'lcn', 'northeastoutside')


roi = 1; % left motor
signal = 2; %
whole_brain_ts =  cell2mat(signals_raw(signal,:));
roi_voxels_ts = whole_brain_ts(I_roi{roi},:);
roi_ts = mean(roi_voxels_ts, 1);
roi_ts_z = zscore(roi_ts);
plot(roi_ts_z);






%
%%%
%%--------------------------------------------------------------------------
%% OFFLINE ANALYSIS/STEPS
%%--------------------------------------------------------------------------
%% Create some comparative plots
%tSNR = mean(kalmanProcTimeSeries, 2)./std(kalmanProcTimeSeries, 0, 2);
%ts = 1:150;
%roi = 1;
%% roi_sigs = [rawTimeSeries(roi,:); glmProcTimeSeries(roi,:); kalmanProcTimeSeries(roi,:); scalProcTimeSeries(roi,:); NFB_disp(roi1,:)];
%% roi_sigs = [glmProcTimeSeries(roi,:); kalmanProcTimeSeries(roi,:); scalProcTimeSeries(roi,:)];
%figure;
%rois = [1,2];
%for roi = 1:numel(rois)
%    roi_sigs = [glmProcTimeSeries(rois(roi),:); kalmanProcTimeSeries(rois(roi),:); scalProcTimeSeries(rois(roi),:)];
%    ax = subplot(3,1,roi)
%    hold on;
%    for i = 1:numel(roi_sigs)
%        plot(ts, roi_sigs(i,:), 'LineWidth', 2)
%    end
%    legend(gca, 'glmProc', 'kalmanProc', 'scaledProc');
%end
%
%%%
%
%figure;
%rois = [1,2,3];
%roi = 1;
%% rawTimeSeries = nan(N_roi, Ndyn);
%% rawTimeSeriesREF = nan(N_roi, Ndyn);
%% glmProcTimeSeries = nan(N_roi, Ndyn);
%% glmPhysProcTimeSeries = nan(N_roi, Ndyn);
%% kalmanProcTimeSeries = nan(N_roi, Ndyn);x
%
%% scalProcTimeSeries = nan(N_roi, Ndyn);
%
%sigs = [convolved_task_design'; rawTimeSeries(rois(roi),:); glmProcTimeSeries(rois(roi),:); kalmanProcTimeSeries(rois(roi),:); scalProcTimeSeries(rois(roi),:); NFB_disp(rois(roi),:)];
%axa = subplot(2,1,1);
%pa = plot(axa, ts, sigs([1 3 4],:), 'LineWidth', 2);
%axb = subplot(2,1,2);
%pb = plot(axb, ts, sigs([1 3 4 5],:), 'LineWidth', 2);
%legend(axa, 'Task', 'glmProc', 'kalmanProc');
%legend(axb, 'Task', 'glmProc', 'kalmanProc', 'scaledProc');
%
%%%
%
%
%ROI_names = {'Left motor','Right motor','GM','WM','CSF','Brain'};
%% ROI_names = cell(1,N_ROIs);
%% ROI_names{1} = 'Left vROI';
%% ROI_names{2} = 'Right vROI';
%% r = 2;
%% ROI_names{r+1} = 'GM';
%% ROI_names{r+2} = 'WM';
%% ROI_names{r+3} = 'CSF';
%% ROI_names{r+4} = 'Brain';
%% ROI_names{r+5} = 'Slice28';
%% ROI_names = cellstr(ROI_names)
%% figure;
%% axa = subplot(1,1,1);
%% hold(axa, 'on')
%% tsdata = rawTimeSeries;
%% % tsdata = glmProcTimeSeries;
%% for roi = 1:numel(ROI_img)
%%     pa = plot(axa, ts, tsdata(roi,:), 'LineWidth', 2);
%%     pa.DisplayName = ROI_names{roi};
%% end
%% legend(axa, ROI_names)
%%
%% rts = tsdata';
%%
%% XX = corrcoef(rts);
%% XXl = tril(XX);
%% figure; imagesc(XXl);
%% ax = gca;
%% ax.YTickLabel = ROI_names;
%% ax.XTickLabel = ROI_names;
%% colorbar
%
%figure;
%axa = subplot(1,1,1);
%hold(axa, 'on')
%% tsdata = rawTimeSeries;
%tsdata = glmProcTimeSeries;
%for roi = 1:numel(ROI_img)
%    pa = plot(axa, ts, tsdata(roi,:), 'LineWidth', 2);
%    pa.DisplayName = ROI_names{roi};
%end
%legend(axa, ROI_names)
%
%rts = tsdata';
%
%XX = corrcoef(rts);
%XXl = tril(XX);
%figure; imagesc(XXl);
%ax = gca;
%ax.YTickLabel = ROI_names;
%ax.XTickLabel = ROI_names;
%colorbar
%
%
%%% Create montage of EPI with ROI(s) as overlay
%% F2D_mean = mean(rF{1}, 2);
%% F3D_mean = reshape(F2D_mean, Nx, Ny, Nz);
%% montage_EPI = onlineBrain_createMontage(F3D_mean, 5, 1, 'Mean EPI (whole image)', 'gray', 'off');
%% f1 = figure;
%% im1 = imagesc(montage_EPI.whole_img);
%% colormap('gray');
%% colorbar;
%% ax = gca;
%% hold(ax, 'on')
%%
%% rois = [1,2,7];
%% roi_color = {'m', 'g', 'r'};
%% [Nimx, Nimy] = size(montage_EPI.whole_img);
%% oo = ones(Nimx, Nimy);
%% zz = zeros(Nimx, Nimy);
%% magenta = cat(3, oo, zz, oo);
%% green = cat(3, zz, oo, zz);
%% red = cat(3, oo, zz, zz);
%% overlay_color = {magenta, green, red};
%% for roi = 1:numel(rois)
%%     montage_roi = onlineBrain_createMontage(ROI_img{rois(roi)}, 5, 1, 'Mask image', 'gray', 'off');
%%     bound_whole_bin = bwboundaries(montage_roi.whole_img);
%%     Nblobs_bin = numel(bound_whole_bin);
%%     for b = 1:Nblobs_bin
%%         p5 = plot(ax, bound_whole_bin{b,1}(:,2), bound_whole_bin{b,1}(:,1), roi_color{roi}, 'LineWidth', 1);
%%     end
%%     imC = imagesc(ax, overlay_color{roi});
%%     set(imC, 'AlphaData', 0.2*montage_roi.whole_img);
%% end
%% hold(ax, 'off');
%
%%% EXTRA CODE - DO NOT DELETE
%% montage_EPI = onlineBrain_createMontageX(F3D_mean, 5, 1, 'Mean EPI (whole image)', 'gray', 'off');
%% f1 = figure;
%% im1 = imagesc(montage_EPI.whole_img);
%
%%%
%% %% Code in progress
%%
%% % Clicking on the brain to select ROI
%%
%% function mouseClick(~,~, bound1, bound2, Nn)
%% disp('kaas')
%% point = get(gca,'CurrentPoint')
%% if (point(1,1) < 0) || (point(1,2) < 0)
%% else
%%     xval = round(point(1,1))
%%     yval = round(point(1,2))
%%     bound1{Nn,1}{1,1}(:,2)
%%     bound1{Nn,1}{1,1}(:,1)
%%
%%     in1 = inpolygon(xval, yval, bound1{Nn,1}{1,1}(:,2), bound1{Nn,1}{1,1}(:,1))
%%     in2 = inpolygon(xval, yval, bound2{Nn,1}{1,1}(:,2), bound2{Nn,1}{1,1}(:,1))
%% end
%%
%% end
%
%% STEP 11: NEUROFEEDBACK SIGNAL CALCULATION / PRESENTATION
%    % OpenNFT does:
%    %   - before starting: initDisplayData % displayData contains fields:
%    %       {'feedbackType', 'condition', 'dispValue', 'Reward',
%    %       'displayStage','displayBlankScreen', 'iteration'};
%    %   - each iteration, from Python: self.displayData = self.eng.nfbCalc(self.iteration, self.displayData, dcmTagLE, dcmOppLE, True, nargout=1)
%    %   - nfbCalc:
%    %       - condition (vectEncCond(indVolNorm)): 1 = index in baseline block; 2 = index in task block
%    %       - get current NF block; and get first index of current NF block
%    %       - get reference baseline in cumulated way across the RUN, or
%    %       other way (get indices)
%    %       - In all ROIs, calculate/get:
%    %           - mBas = median of cumulative baseline;
%    %           - mCond = current value of scalProcTimeSeries;
%    %           - norm_percValues(indRoi) = mCond - mBas;
%    %       - tmp_fbVal = median(norm_percValues);
%    %       - dispValue = round(10000 * tmp_fbVal) /100;
%    %       - keep dispValue between 1 and 100
%
%    %%
%    % Realign (WB)
%% Reslice (WB)
%% (Multi-echo estimation and combination) (WB)
%% Smoothing (WB)
%% iGLM (WB)
%% ROI ==> rawTimeSeries
%% cGLM ==> glmProcTimeSeries
%% Kalman ==> kalmanProcTimeSeries
%% Scaling ==> scalProcTimeSeries
% NFB 