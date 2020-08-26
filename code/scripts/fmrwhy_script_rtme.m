
%--------------------------------------------------------------------------
% DEFINITION
%--------------------------------------------------------------------------
%
%
%
% TODO: NOTE - EVERYTHING IS READ IN AND SAVED WITH SPM_VOL (ETC) AND NOT USING NII_TOOL
%--------------------------------------------------------------------------
% DATA INITIALIZATION FOR ONLINE USE
%--------------------------------------------------------------------------

% Here we mostly create empty variables to be used/populated in real-time
% We also set up all the axes for real-time display
N_roi = N_ROIs;
F = cell(Ne,1);
F_dyn_denoised = F;
rF = F;
srF = nan(Nx*Ny*Nz,Ndyn);
realign_params = F;
fdyn_fn = F;
F_dyn_img = F;
currentVol = F;

I_activeVoxels = cell(1,Ndyn);
NEG_e2n = cell(1,Ndyn);
base3D = cell(Ne,1);

% fMRI data might be multi-echo (i.e. multiple images per volume, not the
% case in the sample data). Create variables per echo.
R = cell(Ne,1);
reslVol = cell(Ne,1);
for e = 1:Ne
    F{e} = nan(Nx*Ny*Nz,Ndyn);
    F_dyn_denoised{e} = nan(Nx*Ny*Nz,Ndyn);
    rF{e}  = nan(Nx*Ny*Nz,Ndyn);
    R{e}(1,1).mat = funcref_spm.mat;
    R{e}(1,1).dim = funcref_spm.dim;
    R{e}(1,1).Vol = funcref_3D;
end
% More Multi-echo variables
if Ne > 1
    e_ref = str2double(options.template_echo);
    X = [ones(Ne,1) -TE(:)];
    base = zeros(Nvox,Ndyn);
    S0 = base;
    T2star = base;
    T2star_corrected = T2star;
    T2star_img = zeros(Nx, Ny, Nz, Ndyn);
    T2star_pv_img = zeros(Nx, Ny, Nz, Ndyn);
    S0_img = zeros(Nx, Ny, Nz, Ndyn);
    S0_pv_img = zeros(Nx, Ny, Nz, Ndyn);
    S_combined_img = zeros(Nx, Ny, Nz, Ndyn);
    S_combined_img2 = zeros(Nx, Ny, Nz, Ndyn);
    S0_pv = base;
    T2star_pv = base;
    T2star_pv_corrected = T2star_pv;
    S0_pv_corrected = S0_pv;
    combined_te_pre = base;
    combined_t2s_pre = base;
    combined_tsnr_pre = base;
    combined_t2s_rt = base;
end

signals_raw_3D = {};
signals_raw = {};
signals_smoothed = {};
initLim = {};
for sig = 1:7
    signals_raw_3D{sig} = nan(Nx, Ny, Nz, Ndyn);
    signals_raw{sig} = nan(Nvox, Ndyn);
    signals_smoothed_3D{sig} = nan(Nx, Ny, Nz, Ndyn);
    signals_smoothed{sig} = nan(Nvox, Ndyn);
    rawTimeSeries{sig} = nan(N_roi, Ndyn);
    displRawTimeSeries{sig} = nan(N_roi, Ndyn);
    rawTimeSeriesREF{sig} = nan(N_roi, Ndyn);
    glmProcTimeSeries{sig} = nan(N_roi, Ndyn);
    glmPhysProcTimeSeries{sig} = nan(N_roi, Ndyn);
    kalmanProcTimeSeries{sig} = nan(N_roi, Ndyn);
    scalProcTimeSeries{sig} = nan(N_roi, Ndyn);
    posMin{sig} = nan(N_roi, Ndyn);
    posMax{sig} = nan(N_roi, Ndyn);
    mposMax{sig} = nan(1, Ndyn);
    mposMin{sig} = nan(1, Ndyn);
    NFB{sig} = nan(N_roi, Ndyn);
    NFB_disp{sig} = nan(N_roi, Ndyn);
    norm_percValues{sig} = nan(N_roi, Ndyn);
end

%%
%--------------------------------------------------------------------------
% SETUP OF VISUALISATIONS
%--------------------------------------------------------------------------

%% Create figure for real-time visualisation
%fig = figure('units','normalized','outerposition',[0 0 1 1]);
%fig_r = 4;
%fig_c = 10;
%fig_tiles = (reshape(1:fig_r*fig_c, fig_c, fig_r))';
%xx = 1:Ndyn;
%
%% 'red' 3D array for tvalue-statistics display
%oo = ones(Nx,Ny);
%zz = zeros(Nx,Ny);
%red = cat(3, oo, zz, zz);
%tmap = zeros(size(funcref_3D));
%tmap_rot = rot90(squeeze(tmap(:,:,Nslice)),rotateVal);
%% Plot background functional image
%funcref_3D_masked = zeros(Nx,Ny,Nz);
%funcref_3D_masked(I_mask) = funcref_3D(I_mask);
%montage_EPI = onlineBrain_createMontage(funcref_3D, 5, 1, 'REF EPI (whole image)', 'gray', 'off');
%% Show:
%%   - Montage, or
%%   - Slice
%rois = [1,2,3];
%roi1 = rois(1);
%roi2 = rois(2);
%ax1 = subplot(fig_r,fig_c,fig_tiles(1:4,1:4));
%showMontage = 0;
%if showMontage
%    % Create montage of refernce EPI with ROI(s) as overlay
%    im1 = imagesc(ax1, montage_EPI.whole_img);
%    colormap('gray');
%    colorbar;
%    hold(ax1, 'on')
%    roi_color = {'m', 'g', 'r'};
%    [Nimx, Nimy] = size(montage_EPI.whole_img);
%    oo = ones(Nimx, Nimy);
%    zz = zeros(Nimx, Nimy);
%    magenta = cat(3, oo, zz, oo);
%    green = cat(3, zz, oo, zz);
%    red = cat(3, oo, zz, zz);
%    overlay_color = {magenta, green, red};
%    for roi = 1:numel(rois)
%        montage_roi = onlineBrain_createMontage(ROI_img{rois(roi)}, 5, 1, 'Mask image', 'gray', 'off');
%        bound_whole_bin = bwboundaries(montage_roi.whole_img);
%        Nblobs_bin = numel(bound_whole_bin);
%        for b = 1:Nblobs_bin
%            p5 = plot(ax1, bound_whole_bin{b,1}(:,2), bound_whole_bin{b,1}(:,1), roi_color{roi}, 'LineWidth', 1);
%        end
%        imC = imagesc(ax1, overlay_color{roi});
%        set(imC, 'AlphaData', 0.2*montage_roi.whole_img);
%    end
%    hold(ax1, 'off');
%else
%    im1 = imagesc(rot90(squeeze(funcref_3D_masked(:,:, floor(Nz/2))),rotateVal)); colormap gray; hold on;
%    im2 = imagesc(red);
%    set(im2, 'AlphaData', tmap_rot);
%    % Plot ROI boundaries as overlay on background
%    Nblobs = numel(bound{roi1}{Nslice,1});
%    for b = 1:Nblobs
%        ax = plot(bound{roi1}{Nslice,1}{b,1}(:,2), bound{roi1}{Nslice,1}{b,1}(:,1), 'm', 'LineWidth', 2);
%        rotate(ax,rotateDir,rotateDeg);
%    end
%    Nblobs = numel(bound{roi2}{Nslice,1});
%    for b = 1:Nblobs
%        ax = plot(bound{roi2}{Nslice,1}{b,1}(:,2), bound{roi2}{Nslice,1}{b,1}(:,1), 'g', 'LineWidth', 2);
%        rotate(ax,rotateDir,rotateDeg);
%    end
%    hold off;
%end
%% Plot time series: task, convolved, task, processed, filtered, scaled
%ax2 = subplot(fig_r,fig_c,fig_tiles(1,5:8)); imagesc(task_design); colormap gray; hold on;
%ln1 = line([1 1], [0.5 1.5], 'Color','red', 'LineWidth', 2); hold off; title('Task/NF design');
%ax3 = subplot(fig_r,fig_c,fig_tiles(2,5:8)); pl1 = plot(1:Ndyn, convolved_task_design, 'LineWidth', 2); xlim([1 Ndyn]); title('HRF convolved task design');
%ax4 = subplot(fig_r,fig_c,fig_tiles(3,5:8)); pl2 = plot(1:Ndyn, rawTimeSeries{1}(1,:), 'LineWidth', 2, 'color', 'm'); xlim([1 Ndyn]); title('Raw NF signal (ROIs: magenta=anat, green=func)'); hold on;
%pl22 = plot(1:Ndyn, rawTimeSeries{1}(2,:), 'LineWidth', 2, 'color', 'g');
%ax5 = subplot(fig_r,fig_c,fig_tiles(4,5:8)); pl3 = plot(1:Ndyn, kalmanProcTimeSeries{1}(1,:), 'LineWidth', 2, 'color', 'm'); xlim([1 Ndyn]); title('Denoised NF signal (ROIs: magenta=anat, green=func)'); hold on;
%pl32 = plot(1:Ndyn, kalmanProcTimeSeries{1}(2,:), 'LineWidth', 2, 'color', 'g');
%set(ax4, 'color','k');
%set(ax5, 'color','k');
%ax6 = subplot(fig_r,fig_c,fig_tiles(1:4,9:10));
%NFB_val = 0; b_handle = bar(NFB_val); ylim([0 100]);
%% Disable callbacks for all objects within the current axes
%set(findall(ax1), 'HitTest', 'off')
%drawnow;

%%
%--------------------------------------------------------------------------
% REAL-TIME ANALYSIS
%--------------------------------------------------------------------------

% Loop through all real-time volumes and conduct iterative analysis
for i = 1:Nt

    tic;
    % Display iteration number
    disp(['i = ' num2str(i)])
    % Skip iteration based on specified amount of initial volumes to skip
    if i <= N_skip
        continue;
    end
    j = i - N_skip;
    
    % STEP 1: LOAD CURRENT VOLUME (all echoes)
    for e = 1:Ne
        % Build dynamic filename
%         fdyn_fn{e} = [func_dir filesep sub '_task-emotion_run-1_echo-2.nii,' num2str(i)];
%         fdyn_fn{e} = [func_dir filesep sub '_task-emotion_run-2_echo-2.nii,' num2str(i)];
        fdyn_fn{e} = fullfile(options.sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_echo-' num2str(e) '_bold_' sprintf('%05d',j) '.nii']);
%         fdyn_fn{e} = [func_dir filesep sub '_task-motor_run-2_echo-2.nii,' num2str(i)];
        currentVol{e} = spm_vol(fdyn_fn{e});
        F_dyn_img{e} = spm_read_vols(currentVol{e}); % Read image volume with SPM. This is the raw/unprocessed image
        F{e}(:,j) = F_dyn_img{e}(:);
    end
    
    % (NOTE: SLICE TIME CORRECTION IS NOT DONE HERE, BUT OFTEN FORMS PART OF A STANDARD FMRI PROCESSING PIPELINE)

    % Current best practice is to realign one echo timeseries and then use
    % the resulting realignment paramaters to calculate and apply a single
    % affine transformation matrix to all echo volumes. Refer to tedana
    % docs

    % Method 1 TODO test vs method 2
    R{e_ref}(2,1).mat = currentVol{e_ref}.mat;
    R{e_ref}(2,1).dim = currentVol{e_ref}.dim;
    R{e_ref}(2,1).Vol = F_dyn_img{e_ref};

    % realign (FROM OPENNFT: preprVol.m)
    [R{e_ref}, A0, x1, x2, x3, wt, deg, b, nrIter] = fmrwhy_realtime_realign(R{e_ref}, flagsSpmRealign, i, N_start, A0, x1, x2, x3, wt, deg, b);
    % get realignment parameters from affine matrices (this does the opposite of spm_matrix, apparently) TODO: check if this does what it is supposed to do
    tmpMCParam = spm_imatrix(R{e_ref}(2,1).mat / R{e_ref}(1,1).mat);
    if (i == N_start)
        offsetMCParam = tmpMCParam(1:6);
    end
    MP(j,:) = tmpMCParam(1:6) - offsetMCParam;

    % Reslice (FROM OPENNFT: preprVol.m) - method 1:
    reslVol{e_ref} = fmrwhy_realtime_reslice(R{e_ref}, flagsSpmReslice, currentVol{e_ref});
    rF{e_ref}(:,j) = reslVol{e_ref}(:);
    for e = 1:Ne
        if e == e_ref
            continue;
        end
%        F{e}(:,:,:,j) = F_dyn_img{e};
        Pm = zeros(12,1);
        Pm(1:6) = MP(j, :);
        orig_mat = currentVol{e}.mat;
        rigid_mat = spm_matrix(Pm, 'T*R');
        trans_mat = rigid_mat * orig_mat;
        R{e}(2,1).dim = currentVol{e}.dim;
        R{e}(2,1).Vol = F_dyn_img{e};
        R{e}(2,1).mat = trans_mat;
        reslVol{e} = fmrwhy_realtime_reslice(R{e}, flagsSpmReslice, currentVol{e});
        rF{e}(:,j) = reslVol{e}(:);
    end
    toc;

%
%    % STEP 2 + 3: REALIGN AND RESLICE TO REFERENCE VOLUME
%    % First realign template echo to template volume
%    % Method 2
%    R(2,1).mat = currentVol{e}.mat;
%    R(2,1).dim = currentVol{e}.dim;
%    R(2,1).Vol = F_dyn_img{e};
%    % realign (FROM OPENNFT: preprVol.m)
%    [R, A0, x1, x2, x3, wt, deg, b, nrIter] = spm_realign_rt(R, flagsSpmRealign, i, N_start, A0, x1, x2, x3, wt, deg, b);
%    % MC params (FROM OPENNFT: preprVol.m)
%    tmpMCParam = spm_imatrix(R(2,1).mat / R(1,1).mat);
%    if (i == N_start)
%        offsetMCParam = tmpMCParam(1:6);
%    end
%    motCorrParam(j,:) = tmpMCParam(1:6)-offsetMCParam; % STEPHAN NOTE: I changed indVolNorm to indVol due to error, not sure if this okay or wrong?
%    MP(j,:) = motCorrParam(j,:);
%    % Reslice (FROM OPENNFT: preprVol.m)
%    reslVol = spm_reslice_rt(R, flagsSpmReslice);
%    rF{e}(:,j) = reslVol(:);
%
%    for e = 1:Ne
%        if e == str2double(options.template_echo)
%            continue;
%        end
%
%
%    end


    % STEP 4: MULTI-ECHO PARAMETER ESTIMATION AND COMBINATION (not for sample data)
    if Ne > 1
        % If option selected to use all echoes (i.e. use_echo = 0) continue
        % with multi-echo parameter estimation and combination, else use
        % specified echo (e.g use_echo = 2)
        if use_echo == 0

            me_params = fmrwhy_realtime_estimateMEparams(reslVol, options.TE, I_mask);
            T2star_pv(:,j) = reshape(me_params.T2star_3D, Nx*Ny*Nz, 1);
            S0_pv(:,j) = reshape(me_params.S0_3D, Nx*Ny*Nz, 1);
            T2star_pv_corrected(:,j) = reshape(me_params.T2star_3D_thresholded, Nx*Ny*Nz, 1);
            S0_pv_corrected(:,j) = reshape(me_params.S0_3D_thresholded, Nx*Ny*Nz, 1);

            func_data = zeros(Nx, Ny, Nz, Ne);
            for e = 1:Ne
                func_data(:,:,:,e) = reslVol{e};
            end
            combined_t2s_pre_3D = fmrwhy_me_combineEchoes(func_data, options.TE, 0, 1, t2star_img);
            combined_t2s_rt_3D = fmrwhy_me_combineEchoes(func_data, options.TE, 0, 1, me_params.T2star_3D);
            combined_tsnr_pre_3D = fmrwhy_me_combineEchoes(func_data, options.TE, 0, 2, tsnr_data);
            combined_te_pre_3D = fmrwhy_me_combineEchoes(func_data, options.TE, 0, 3, options.TE);
            combined_t2s_pre(:,j) = combined_t2s_pre_3D(:);
            combined_t2s_rt(:,j) = combined_t2s_rt_3D(:);
            combined_tsnr_pre(:,j) = combined_tsnr_pre_3D(:);
            combined_te_pre(:,j) = combined_te_pre_3D(:);

            signals_raw_3D{1}(:,:,:,j) = reslVol{2};                       % Echo 2
            signals_raw_3D{2}(:,:,:,j) = combined_tsnr_pre_3D;             % Pre-tSNR combined
            signals_raw_3D{3}(:,:,:,j) = combined_t2s_pre_3D;              % Pre-T2* combined
            signals_raw_3D{4}(:,:,:,j) = combined_te_pre_3D;               % Pre-TE combined
            signals_raw_3D{5}(:,:,:,j) = combined_t2s_rt_3D;               % RT-T2* combined
            signals_raw_3D{6}(:,:,:,j) = me_params.T2star_3D_thresholded;  % T2* FIT
            signals_raw_3D{7}(:,:,:,j) = me_params.S0_3D_thresholded;      % S0 FIT

            signals_raw{1}(:,j) = reslVol{2}(:);
            signals_raw{2}(:,j) = combined_tsnr_pre(:,j);
            signals_raw{3}(:,j) = combined_t2s_pre(:,j);
            signals_raw{4}(:,j) = combined_te_pre(:,j);
            signals_raw{5}(:,j) = combined_t2s_rt(:,j);
            signals_raw{6}(:,j) = T2star_pv_corrected(:,j);
            signals_raw{7}(:,j) = S0_pv_corrected(:,j);

            rf = combined_tsnr_pre_3D;
        else
            rf = rF{use_echo}(:,j);
        end
    else
        % if single-echo, use first volume in rF cell array
        rf = rF{1}(:,j);
    end

    % STEP 5: SMOOTH REALIGNED VOLUME
    % Using OpenNFT functionality and SPM
    for sig = 1:numel(signals_raw_3D)
        srf = zeros(Nx, Ny, Nz);
        gKernel = smoothing_kernel ./ dicomInfoVox;
        spm_smooth(squeeze(signals_raw_3D{sig}(:,:,:,j)), srf, gKernel);
        signals_smoothed_3D{sig}(:,:,:,j) = srf;
        signals_smoothed{sig}(:,j) = srf(:);
        srF(:,j) = srf(:);
    end

%    % STEP 6: AR(1) FILTERING OF SMOOTHED VOLUME (for neufep: false)
%    if iglmAR1
%        if j == 1
%            % initalize first AR(1) volume
%            asrF(:,j) = (1 - aAR1) * srF(:,j);
%        else
%            asrF(:,j) = srF(:,j) - aAR1 * asrF(:,j-1);
%        end
%    else
%        asrF(:,j) = srF(:,j);
%    end
%
%    % STEP 7: iGLM FOR VOLUME (for neufep: false)
%    if isIGLM
%        % Scaling settings
%        if fLockedTempl
%            if j == 1
%                % Only set the scaling settings based on first iteration
%                max_smReslVol = max(asrF(:,j));
%                min_smReslVol = min(asrF(:,j));
%                normSmReslVol = (asrF(:,j)-min_smReslVol) / (max_smReslVol-min_smReslVol);
%            end
%        else
%            % Update scaling settings on each iteration
%            max_smReslVol = max(asrF(:,j));
%            min_smReslVol = min(asrF(:,j));
%            normSmReslVol = (asrF(:,j)-min_smReslVol) / (max_smReslVol-min_smReslVol);
%        end
%
%        % Constant regressor is always included
%        constRegr = constRegrFull(1:j);
%        if isRegrIGLM
%            % Create empty nuisance regressors
%            motRegr = [];
%            linRegr = [];
%            highPassRegr = [];
%            % Set regressor content if they need to be included in design
%            % matrix
%            if isMotionRegr
%                if Ne > 1
%                    motRegr = zscore(MP(1:j,:));
%                else
%                    motRegr = zscore(MP(1:j,:));
%                end
%            end
%            if isLinRegr
%                linRegr = linRegrFull(1:j);
%            end
%            if isHighPass
%                highPassRegr = cosine_basis_set(1:j, :);
%            end
%            % Construct design matrix without task/baseline conditions, i.e.
%            % including nuisance and constant regressors
%            tmpRegr = horzcat(motRegr, linRegr, highPassRegr, constRegr);
%        else
%            tmpRegr = constRegr;
%        end
%        nrBasFctRegr = size(tmpRegr, 2);
%
%        % AR(1) for regressors of no interest
%        if iglmAR1
%            tmpRegr = arRegr_opennft(aAR1,tmpRegr);
%        end
%        % combine with prepared basFct design regressors
%        basFctRegr = [basFct(1:j,:), tmpRegr];
%
%        % estimate iGLM
%        [idxActVoxIGLM, dyntTh, tTh, Cn, Dn, s2n, tn, neg_e2n] = ...
%            onlineBrain_iGLM(Cn, Dn, s2n, tn, asrF(:,j), j, ...
%            (nrBasFct+nrBasFctRegr), tContr, basFctRegr, pVal, ...
%            dyntTh, tTh, spmMaskTh);
%
%        % catch negative iGLM estimation error message for log
%        NEG_e2n{j} = neg_e2n;
%        if ~isempty(neg_e2n)
%            disp('HERE THE NEGATIVE e2n!!!')
%        end
%        I_activeVoxels{j} = idxActVoxIGLM;
%    else
%        idxActVoxIGLM = [];
%    end
%    % handle empty activation map and division by 0
%    if ~isempty(idxActVoxIGLM) && max(tn) > 0
%        maskedStatMapVect = tn(idxActVoxIGLM);
%        maxTval = max(maskedStatMapVect);
%        statMapVect = maskedStatMapVect;
%        statMap3D(idxActVoxIGLM) = statMapVect;
%        statMap4D{j} = statMap3D;
%    end

    % STEP 8 + 9 + 10: cGLM NUISANCE REGRESSION, KALMAN FILTERING AND SCALING OF SIGNAL IN ROI(s)
%    if isPhysRegr
%        rawTimeSeriesREF(N_ROI_REF,j) = mean(asrF(I_roi{N_ROI_REF},j));
%    end

    for roi = 1:numel(ROI_img)

        for sig = 1:numel(signals_raw_3D)
            rawTimeSeries{sig}(roi,j) = mean(signals_smoothed{sig}(I_roi{roi},j));

            % Limits for scaling
            initLim{sig}(roi) = 0.005*mean(rawTimeSeries{sig}(roi,1:j));

            % Raw for Display (STEPHAN: FIGURE OUT WHY THIS IS DONE)
            displRawTimeSeries{sig}(roi,j) = rawTimeSeries{sig}(roi, j)-rawTimeSeries{sig}(roi, 1);

            % To avoid NaNs given algnment to zero, see preprVol()
            motCorrParam(1,:) = 0.00001;

            % Get full time series up to current iteration
            tmp_rawTimeSeries = rawTimeSeries{sig}(roi, 1:j)';
            % tmp_rawTimeSeriesREF = rawTimeSeriesREF(roi, 1:j)';

            % Time-series AR(1) filtering
            if cglmAR1
                % initalize first AR(1) value
                if j == 1
                    tmp_rawTimeSeriesAR1(roi,j) = (1 - aAR1) * tmp_rawTimeSeries(j);
                else
                    tmp_rawTimeSeriesAR1(roi,j) = tmp_rawTimeSeries(j) - aAR1 * tmp_rawTimeSeriesAR1(roi,j-1);
                end
                % replace raw ime-series with AR(1) time-series
                clear tmp_rawTimeSeries
                tmp_rawTimeSeries = tmp_rawTimeSeriesAR1(roi, :)';
            end

            % Setup design matrix regressors
            constRegrC = constRegrFull(1:j); % Constant regressor is always included
            linRegrC = [];
            motRegrC = [];
            designRegrC = [];
            physRegrC = [];

            % Step-wise addition of regressors, step = total nr of regressors,
            % which may require a justification for particular project
            regrStep = nrRegrToCorrect + nrRegrDesign; % 1 + 1 + 6 + 1 = 9
            if j < regrStep
                % only include constant regressor
            elseif (j >= regrStep) && (j < 2*regrStep)
                % include constant and linear regressors
                linRegrC = linRegrFull(1:j);
            else %(j >= 2*regrStep)
                % include constant, linear and motion correction regressors
                linRegrC = linRegrFull(1:j);
                if Ne > 1
                    motRegrC = zscore(MP(1:j,:));
                else
                    motRegrC = zscore(MP(1:j,:));
                end
                if isPhysRegr
                    physRegrC = rawTimeSeriesREF(N_ROI_REF,1:j)';
                end
            end

            % Concatenate regressors of no interest
            tmpRegr = horzcat(constRegrC, linRegrC, motRegrC);
            if isPhysRegr
                tmpRegr = horzcat(constRegrC, linRegrC, motRegrC, physRegrC);
            end

            % AR(1) for regressors of no interest
            if cglmAR1
                tmpRegr = arRegr_opennft(aAR1,tmpRegr); %TODO, rename and move this function?
            end

            % Create final design matrix and estimate GLM parameters
            if j < 3*regrStep
                % estimate GLM parameters for case where task regressor is not
                % included
                cX = tmpRegr;
                beta = pinv(cX) * tmp_rawTimeSeries;
                tmp_glmProcTimeSeries = (tmp_rawTimeSeries - cX * beta)';
            else
                % First include task regressor into design matrix
                cX = [tmpRegr spmDesign(1:j,:)];
                beta = pinv(cX) * tmp_rawTimeSeries;
                tmp_glmProcTimeSeries = (tmp_rawTimeSeries - cX * [beta(1:end-1); zeros(1,1)])';
            end

            glmProcTimeSeries{sig}(roi, j) = tmp_glmProcTimeSeries(end);

            % Modified Kalman low-pass filter + spike identification & correction
            tmpStd = std(glmProcTimeSeries{sig}(roi,1:j));
            S(roi).Q = tmpStd^2;
            S(roi).R = 1.95*tmpStd^2;
            kalmThreshold = 0.9*tmpStd;
            [kalmanProcTimeSeries{sig}(roi,j), S(roi), fPositDerivSpike(roi), fNegatDerivSpike(roi)] = ...
                onlineBrain_modifKalman(kalmThreshold, glmProcTimeSeries{sig}(roi,j), S(roi), fPositDerivSpike(roi), fNegatDerivSpike(roi));

            % Scaling: TODO: decide if need to SKIP KALMAN FILTERING ==> PUT GLMPROCTIMESERIES INTO SCALING FUNCTION BELOW
            slWind = basBlockLength * nrBlocksInSlidingWindow;
            [scalProcTimeSeries{sig}(roi, j), tmp_posMin{sig}(roi), tmp_posMax{sig}(roi)] = ...
                onlineBrain_scaleTimeSeries(kalmanProcTimeSeries{sig}(roi,1:j), j, slWind, basBlockLength, initLim{sig}(roi), vectEncCond(1:j), tmp_posMin{sig}(roi), tmp_posMax{sig}(roi));
            posMin{sig}(roi,j)=tmp_posMin{sig}(roi);
            posMax{sig}(roi,j)=tmp_posMax{sig}(roi);

        end
    end

    for sig = 1:numel(signals_raw_3D)
        % calcualte average limits for 2 ROIs, e.g. for bilateral NF
        % NF extensions with >2 ROIs requires an additional justification
        mposMax{sig}(j)= mean(posMax{sig}(:, j));
        mposMin{sig}(j)= mean(posMin{sig}(:, j));

        % STEP 11: NEUROFEEDBACK SIGNAL CALCULATION / PRESENTATION
        if baseline_design(j) == 1
            % If the current iteration is in a baseline block, feedback value
            % is zero for all ROIs
            NFB{sig}(:,j) = 0;
            NFB_disp{sig}(:,j) = 0;
        else
            % If the current iteration is in a task block, feedback value
            % is calculated as PSC of current value compared to cumulative
            % basline mean/median. This is done per ROI
            i_bas = I_baseline(I_baseline<j);
            for roi = 1:numel(ROI_img)
                mBas = median(scalProcTimeSeries{sig}(roi,i_bas));
                mCond = scalProcTimeSeries{sig}(roi,j);
                norm_percValues{sig}(roi, j) = mCond - mBas;
                % tmp_fbVal = median(norm_percValues); ==> OpenNFT calculates
                % median over ROIs, e.g. when interested in signal in multiple
                % ROIs like bilateral occipital cortices
                NFB{sig}(roi,j) = norm_percValues{sig}(roi, j);
                NFB_disp{sig}(roi,j) = round(10000 * NFB{sig}(roi,j)) /100;
                % [1...100], for Display
                if NFB_disp{sig}(roi,j) < 1
                    NFB_disp{sig}(roi,j) = 1;
                end
                if NFB_disp{sig}(roi,j) > 100
                    NFB_disp{sig}(roi,j) = 100;
                end
            end
        end
    end

    % STEP 12: UPDATE PLOTS
    % TODO: ROIPOLY to draw a polygon for user-specified roi (e.g. extract signal)
    % TODO: INPOLYGON to check if buttonpress is inside mask/ROI boundary
%    if showMontage
%    else
%        if ~isempty(idxActVoxIGLM)
%            tmap = statMap4D{j};
%            tmap_masked = zeros(size(tmap));
%            tmap_masked(I_mask) = tmap(I_mask);
%            tmap_rot = rot90(squeeze(tmap_masked(:,:,Nslice)),rotateVal);
%            normA = tmap_rot - min(tmap_rot(:));
%            normA = normA ./ max(normA(:));
%            set(im2, 'AlphaData', normA);
%        end
%    end
%    set(ln1, 'XData', [j j]);
%    plt_sig = 2;
%    set(pl2, 'YData', rawTimeSeries{plt_sig}(roi1,:));
%    set(pl22, 'YData', rawTimeSeries{plt_sig}(roi2,:));
%    set(pl3, 'YData', kalmanProcTimeSeries{plt_sig}(roi1,:)); % TODO: change between displaying kalman filtered time series and glm processed time series
%    set(pl32, 'YData', kalmanProcTimeSeries{plt_sig}(roi2,:));
%    set(b_handle, 'Ydata', NFB_disp{plt_sig}(roi2,j))
%    drawnow;
end
