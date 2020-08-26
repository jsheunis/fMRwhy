%--------
%--------
% Script to post-process the real-time fMRI task data
%--------
%--------

tSNR = {};
tSNR_fns = {};
tSNR_MNI_fns = {};
mean_fns = {};
stddev_fns = {};
sig_desc = {'RTecho2', 'RTcombinedTSNR', 'RTcombinedT2STAR', 'RTcombinedTE', 'RTcombinedRTt2star', 'RTt2starFIT', 'RTs0FIT'};
template_fn = options.template_fn;
sub_dir_rt = options.sub_dir_rt;


%--------
% STEP 1: Save functional timeseries to nifti
%--------

rfunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rpreproc_bold.nii']);
for sig = 1:numel(signals_raw_3D)
    new_spm_raw = spm_vol(rfunctional_fn);
    new_spm_smoothed = spm_vol(rfunctional_fn);
    realigned_time_series_fn = fullfile(sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-r' sig_desc{sig} '_bold.nii']);
    smoothed_time_series_fn = fullfile(sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-sr' sig_desc{sig} '_bold.nii']);

    if ~exist(realigned_time_series_fn, 'file')
        for i = 1:numel(new_spm_raw)
            new_spm_raw(i).fname = realigned_time_series_fn;
            new_spm_raw(i).private.dat.fname = realigned_time_series_fn;
            new_spm_smoothed(i).fname = smoothed_time_series_fn;
            new_spm_smoothed(i).private.dat.fname = smoothed_time_series_fn;
            if sig == 6
                new_spm_raw(i).pinfo(1) = 1;
                new_spm_smoothed(i).pinfo(1) = 1;
            end
            spm_write_vol(new_spm_raw(i), signals_raw_3D{sig}(:,:,:,i));
            spm_write_vol(new_spm_smoothed(i), signals_smoothed_3D{sig}(:,:,:,i));
        end
    end
end


%--------
% STEP 2: Save NFB roi timeseries to tsv
%--------

roi_desc = {'anat', 'funcAnat', 'GM', 'WM', 'CSF', 'brain'};
nfb_sigs = {'raw', 'rawDisp', 'glm', 'kalm', 'scal', 'normPerc', 'nfb', 'nfbDisp'}

for roi = 1:numel(roi_desc)
    col_names = [];
    txt_fn = fullfile(sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-' roi_desc{roi} '_ROIsignals.txt']);
    data = [];
    for sig = 1:numel(sig_desc)
        names = {};
        for n = 1:numel(nfb_sigs-1) % exclude S0 time series
            names{n} = [nfb_sigs{n} '_' sig_desc{sig}];
        end
        col_names = [col_names, names];

        data = [data rawTimeSeries{sig}(roi,:)' displRawTimeSeries{sig}(roi,:)' glmProcTimeSeries{sig}(roi,:)' kalmanProcTimeSeries{sig}(roi,:)' scalProcTimeSeries{sig}(roi,:)' norm_percValues{sig}(roi,:)' NFB{sig}(roi,:)' NFB_disp{sig}(roi,:)'];
    end

    data_table = array2table(data,'VariableNames', col_names);
    writetable(data_table, txt_fn, 'Delimiter','\t');
    [status, msg, msgID] = movefile(temp_txt_fn, tsv_fn);

%    dlmwrite(txt_fn, data, 'delimiter', '\t', 'precision', '%1.7e')
%    tsv_fn = fmrwhy_util_saveAsTSV(txt_fn, col_names);
end


%--------
% STEP 3: Calculate mean, stddev, tSNR
%--------

for sig = 1:numel(sig_desc)

    data_4D = signals_raw_3D{sig}; % [Ni x Ny x Nz x Ndyn]
    [Nx, Ny, Nz, Ndyn] = size(data_4D);
    data_2D = reshape(data_4D, Nx*Ny*Nz, Ndyn); %[voxels, time]
    data_2D = data_2D'; %[time, voxels]

    % Remove linear and quadratic trends from data, per voxel
    data_2D_detrended = fmrwhy_util_detrend(data_2D, 2); %[time, voxels]

    % Calculate mean (ignore nan values)
    data_2D_mean = nanmean(data_2D_detrended); %[1, voxels]
    data_3D_mean = reshape(data_2D_mean, Nx, Ny, Nz);

    % Calculate standard deviation
    data_2D_stddev = std(data_2D_detrended); %[1, voxels]
    data_3D_stddev = reshape(data_2D_stddev, Nx, Ny, Nz);

    % Calculate tSNR
    data_2D_tsnr = data_2D_mean./data_2D_stddev; %[1, voxels]
    data_2D_tsnr(data_2D_tsnr<0) = 0; % TODO: double check if this thresholding is fine to do, or give adequate reasoning
    data_3D_tsnr = reshape(data_2D_tsnr, Nx, Ny, Nz);

    % Mask if required
    data_2D_tsnr_masked = zeros(size(data_2D_tsnr)); %[1, voxels]
    data_2D_tsnr_masked(I_mask) = data_2D_tsnr(I_mask); %[1, voxels]
    data_3D_tsnr_masked = reshape(data_2D_tsnr_masked, Nx, Ny, Nz);
    tSNR{sig} = data_3D_tsnr_masked;

    % Save to file, if required
    tSNR_fns{sig} = fullfile(sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-' sig_desc{sig} '_tsnr.nii']);
    tSNR_MNI_fns{sig} = fullfile(sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_space-MNI152_desc-' sig_desc{sig} '_tsnr.nii']);
    mean_fns{sig} = fullfile(sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-' sig_desc{sig} '_mean.nii']);
    stddev_fns{sig} = fullfile(sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-' sig_desc{sig} '_stddev.nii']);
    if ~exist(tSNR_fns{sig})
        fmrwhy_util_saveNiftiOld(tSNR_fns{sig}, data_3D_tsnr, template_fn, sig_desc{sig}, 1)
    end
    if ~exist(mean_fns{sig})
        fmrwhy_util_saveNiftiOld(mean_fns{sig}, data_3D_mean, template_fn, sig_desc{sig}, 1)
    end
    if ~exist(stddev_fns{sig})
        fmrwhy_util_saveNiftiOld(stddev_fns{sig}, data_3D_stddev, template_fn, sig_desc{sig}, 1)
    end

end


%--------
% STEP 4: Compare tSNR maps (perc signal change, distributions, etc)
%--------

% Mask details
masks = fmrwhy_util_loadMasks(bids_dir, sub);
mask_fn = masks.brain_mask_fn;
% tSNR and other filenames
tsnr_pngs = {};
percdiff_pngs = {};
for i = 1:numel(tSNR_fns)
    tsnr_pngs{i} = strrep(tSNR_fns{i}, '.nii', '.png');
    if i > 1
        percdiff_pngs{i-1} = strrep(tSNR_fns{i}, '_tsnr.nii', '_percdiff.png');
    end
end
distr_png = fullfile(sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-tsnrPercdiffRainclouds.png']);
% ROIs
roi_fns = {};
roi_fns{1} = roi_anat_fn;
roi_fns{2} = roi_funcAnat_fn;
if strcmp(task, 'motor')
    compare_roi_txt = {': anatomical ROI - leftMotor', ': functional ROI - leftMotor'};
else
    compare_roi_txt = {': anatomical ROI - bilatAmygdala', ': functional ROI - bilatAmygdala'};
end
%fmrwhy_util_compareTSNR(tsnr_fns, mask_fn, roi_fns, compare_roi_txt, tsnr_saveAs_fns, perc_saveAs_fns, distr_saveAs_fn)
fmrwhy_util_compareTSNRrt(tSNR_fns, mask_fn, roi_fns, compare_roi_txt, tsnr_pngs, percdiff_pngs, distr_png)


% -------
% STEP 5: tSNR image warping
% -------

toTransform_fns = tSNR_fns;
saveAs_transform_fns = tSNR_MNI_fns;
transformation_fn = options.indiv_to_mni_fn;
if ~exist(toTransform_fns{1}, 'file')
    fmrwhy_batch_normaliseWrite(toTransform_fns, transformation_fn, template_fn, saveAs_transform_fns)
end




% -------
% STEP 6: Delineate tSNR values per tissue type and ROI
% -------
roi_desc_txt = {'lmotor', 'bamygdala'};
masks_oriented = fmrwhy_util_loadOrientMasks(bids_dir, sub);
mask_img_oriented = masks_oriented.brain_mask_3D;
%  TSVs with tsnr values extracted per tissue mask (GM, WM, CSF, whole brain) and ROI
for i = 1:numel(tSNR_fns)
    [p_tsnr, frm, rg, dim] = fmrwhy_util_readOrientNifti(tSNR_fns{i});
    tsnr_img = p_tsnr.nii.img(:);
    for j = 1:4
        vals = tsnr_img(masks_oriented.([masks_oriented.field_names{j} '_mask_I']));
        tsnr_output_fn = strrep(tSNR_fns{i}, '_tsnr.nii', ['_' masks_oriented.field_names{j} 'tsnr.tsv']);
        temp_txt_fn = strrep(tsnr_output_fn, '.tsv', '_temp.txt');
        data_table = array2table(vals,'VariableNames', {'tsnr'});
        writetable(data_table, temp_txt_fn, 'Delimiter','\t');
        [status, msg, msgID] = movefile(temp_txt_fn, tsnr_output_fn);
    end

    for k = 1:numel(roi_fns)
        [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(roi_fns{k});
        roi_img = fmrwhy_util_createBinaryImg(p.nii.img, 0.1);
        roi_img_2D = roi_img(:);
        I_roi = find(roi_img_2D);
        masked_tsnr_img = fmrwhy_util_maskImage(p_tsnr.nii.img, roi_img);

        for j = 1:4
            overlap = masks_oriented.([masks_oriented.field_names{j} '_mask_2D']) & roi_img_2D;
            vals = tsnr_img(find(overlap));
            tsnr_output_fn = strrep(tSNR_fns{i}, '_tsnr.nii', ['_' roi_desc_txt{k} masks_oriented.field_names{j} 'tsnr.tsv']);
            temp_txt_fn = strrep(tsnr_output_fn, '.tsv', '_temp.txt');
            data_table = array2table(vals,'VariableNames', {'tsnr'});
            writetable(data_table, temp_txt_fn, 'Delimiter','\t');
            [status, msg, msgID] = movefile(temp_txt_fn, tsnr_output_fn);
        end
    end
end


% -------
% STEP 7: ROI Timeseries plot outputs
% -------

for sig = 1:numel(sig_desc)
    functional_fn = fullfile(sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-sr' sig_desc{sig} '_bold.nii']);

    task_info.TR = options.firstlevel.(task).(['run' run]).sess_params.timing_RT;
    task_info.onsets = options.firstlevel.(task).(['run' run]).plot_params.cond_onset;
    task_info.durations = options.firstlevel.(task).(['run' run]).plot_params.cond_duration;
    task_info.precision = 1;
    trace_info = [];

    roi_fn = roi_anat_fn;
    trace_info = NFB{sig}(1,:);
    saveAs_fn = strrep(functional_fn, '_bold.nii', '_anatROItsplot.png');
    fmrwhy_util_thePlotROIrt(functional_fn, mask_fn, roi_fn, task_info, trace_info, saveAs_fn)

    roi_fn = roi_funcAnat_fn;
    trace_info = NFB{sig}(2,:);
    saveAs_fn = strrep(functional_fn, '_bold.nii', '_funcAnatROItsplot.png');
    fmrwhy_util_thePlotROIrt(functional_fn, mask_fn, roi_fn, task_info, trace_info, saveAs_fn)
end


% TODO:
% - tCNR





%
%figure;
%hold on;
%leg = {};
%[cb] = cbrewer('qual','Set3',12,'pchip');
%colorss = {cb(1,:), cb(5,:), cb(4,:), cb(6,:), cb(8,:), cb(10,:)}
%pp = {};
%for i = 1:6
%    pp{i} = plot(1:Ndyn, kalmanProcTimeSeries{i}(roi2,:), 'color', colorss{i}, 'LineWidth', 1.5)
%end
%
%legend('Echo', 'T2{\ast}-combined', 'tSNR-combined', 'TE-combined', 'rtT2{\ast}-combined', 'rtT2{\ast}')
%




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