%--------
%--------
% Script to post-explore the real-time fMRI task data
%--------
%--------

% -------
% STEP 1: Define run details
% -------
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

tasks = {'motor', 'emotion'};
tasks = {'motor'};
runs = {'1', '2'};
runs = {'1'};
echoes = {'2', 'combinedMEtsnr', 'combinedMEt2star', 'combinedMEte'};
subs = {'001', '002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
subs = {'001'};
ses = '';
echo = echoes{1};

sig_desc = {'RTecho2', 'RTcombinedTSNR', 'RTcombinedT2STAR', 'RTcombinedTE', 'RTcombinedRTt2star', 'RTt2starFIT', 'RTs0FIT'};

for s = 1:numel(subs)

    clearvars -except bids_dir tasks runs echoes subs ses echo s

    sub = subs{s};

    options = struct;
    % Setup fmrwhy BIDS-derivatuve directories on workflow level
    options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);
    % Grab parameters from workflow settings file
    options = fmrwhy_settings_preprocQC(bids_dir, options);
    % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);
    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_defaults_subAnat(bids_dir, sub, options);
    % Multi-echo derivatives
    options.me_dir = fullfile(options.deriv_dir, 'fmrwhy-multiecho');
    options.sub_dir_me = fullfile(options.me_dir, ['sub-' sub]);
    % Create derivatives directory for rt output
    options.rt_dir = fullfile(options.deriv_dir, 'fmrwhy-rt');
    if ~exist(options.rt_dir, 'dir')
        mkdir(options.rt_dir);
    end
    % Create sub directory for rt output
    options.sub_dir_rt = fullfile(options.rt_dir, ['sub-' sub]);
    if ~exist(options.sub_dir_rt, 'dir')
        mkdir(options.sub_dir_rt);
    end

    for t = 1:numel(tasks)
        task = tasks{t};

        for r = 1:numel(runs)
            run = runs{r};

            disp('-----------------------------------------------------')
            disp(['Running RT analysis for: sub-' sub '_task-' task '_run-' run])
            disp('-----------------------------------------------------')

            % Update workflow params with subject functional derivative filenames
            options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);
            realigned_t2star = fullfile(options.sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-rRTt2starFIT_bold.nii']);
            smoothed_t2star = fullfile(options.sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-srRTt2starFIT_bold.nii']);
            realigned_s0 = fullfile(options.sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-rRTs0FIT_bold.nii']);
            smoothed_s0 = fullfile(options.sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-srRTs0FIT_bold.nii']);

            t2star_nii = nii_tool('load', realigned_t2star);
            t2star_4D = t2star_nii.img;
            [Nx, Ny, Nz, Nt] = size(t2star_4D);
            t2star_2D = reshape(t2star_4D, Nx*Ny*Nz, Nt);

            s0_nii = nii_tool('load', realigned_s0);
            s0_4D = s0_nii.img;
            s0_2D = reshape(s0_4D, Nx*Ny*Nz, Nt);


%            for sig = 1:numel(sig_desc)
%                realigned_time_series_fn = fullfile(sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-r' sig_desc{sig} '_bold.nii']);
%                smoothed_time_series_fn = fullfile(sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-sr' sig_desc{sig} '_bold.nii']);
%
%            end

        end
    end
end
%
%
%
%%--------
%% STEP 4: Compare tSNR maps (perc signal change, distributions, etc)
%%--------
%
%% Mask details
%masks = fmrwhy_util_loadMasks(bids_dir, sub);
%mask_fn = masks.brain_mask_fn;
%% tSNR and other filenames
%tsnr_pngs = {};
%percdiff_pngs = {};
%for i = 1:numel(tSNR_fns)
%    tsnr_pngs{i} = strrep(tSNR_fns{i}, '.nii', '.png');
%    if i > 1
%        percdiff_pngs{i-1} = strrep(tSNR_fns{i}, '_tsnr.nii', '_percdiff.png');
%    end
%end
%distr_png = fullfile(sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-tsnrPercdiffRainclouds.png']);
%% ROIs
%roi_fns = {};
%roi_fns{1} = roi_anat_fn;
%roi_fns{2} = roi_funcAnat_fn;
%if strcmp(task, 'motor')
%    compare_roi_txt = {': anatomical ROI - leftMotor', ': functional ROI - leftMotor'};
%else
%    compare_roi_txt = {': anatomical ROI - bilatAmygdala', ': functional ROI - bilatAmygdala'};
%end
%%fmrwhy_util_compareTSNR(tsnr_fns, mask_fn, roi_fns, compare_roi_txt, tsnr_saveAs_fns, perc_saveAs_fns, distr_saveAs_fn)
%fmrwhy_util_compareTSNRrt(tSNR_fns, mask_fn, roi_fns, compare_roi_txt, tsnr_pngs, percdiff_pngs, distr_png)
%
%
%% -------
%% STEP 5: tSNR image warping
%% -------
%
%toTransform_fns = tSNR_fns;
%saveAs_transform_fns = tSNR_MNI_fns;
%transformation_fn = options.indiv_to_mni_fn;
%if ~exist(toTransform_fns{1}, 'file')
%    fmrwhy_batch_normaliseWrite(toTransform_fns, transformation_fn, template_fn, saveAs_transform_fns)
%end
%
%
%
%
%% -------
%% STEP 6: Delineate tSNR values per tissue type and ROI
%% -------
%roi_desc_txt = {'lmotor', 'bamygdala'};
%masks_oriented = fmrwhy_util_loadOrientMasks(bids_dir, sub);
%mask_img_oriented = masks_oriented.brain_mask_3D;
%%  TSVs with tsnr values extracted per tissue mask (GM, WM, CSF, whole brain) and ROI
%for i = 1:numel(tSNR_fns)
%    [p_tsnr, frm, rg, dim] = fmrwhy_util_readOrientNifti(tSNR_fns{i});
%    tsnr_img = p_tsnr.nii.img(:);
%    for j = 1:4
%        vals = tsnr_img(masks_oriented.([masks_oriented.field_names{j} '_mask_I']));
%        tsnr_output_fn = strrep(tSNR_fns{i}, '_tsnr.nii', ['_' masks_oriented.field_names{j} 'tsnr.tsv']);
%        temp_txt_fn = strrep(tsnr_output_fn, '.tsv', '_temp.txt');
%        data_table = array2table(vals,'VariableNames', {'tsnr'});
%        writetable(data_table, temp_txt_fn, 'Delimiter','\t');
%        [status, msg, msgID] = movefile(temp_txt_fn, tsnr_output_fn);
%    end
%
%    for k = 1:numel(roi_fns)
%        [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(roi_fns{k});
%        roi_img = fmrwhy_util_createBinaryImg(p.nii.img, 0.1);
%        roi_img_2D = roi_img(:);
%        I_roi = find(roi_img_2D);
%        masked_tsnr_img = fmrwhy_util_maskImage(p_tsnr.nii.img, roi_img);
%
%        for j = 1:4
%            overlap = masks_oriented.([masks_oriented.field_names{j} '_mask_2D']) & roi_img_2D;
%            vals = tsnr_img(find(overlap));
%            tsnr_output_fn = strrep(tSNR_fns{i}, '_tsnr.nii', ['_' roi_desc_txt{k} masks_oriented.field_names{j} 'tsnr.tsv']);
%            temp_txt_fn = strrep(tsnr_output_fn, '.tsv', '_temp.txt');
%            data_table = array2table(vals,'VariableNames', {'tsnr'});
%            writetable(data_table, temp_txt_fn, 'Delimiter','\t');
%            [status, msg, msgID] = movefile(temp_txt_fn, tsnr_output_fn);
%        end
%    end
%end
%
%
%% -------
%% STEP 7: ROI Timeseries plot outputs
%% -------
%
%for sig = 1:numel(sig_desc)
%    functional_fn = fullfile(sub_dir_rt, ['sub-' sub '_task-' task '_run-' run '_desc-sr' sig_desc{sig} '_bold.nii']);
%
%    task_info.TR = options.firstlevel.(task).(['run' run]).sess_params.timing_RT;
%    task_info.onsets = options.firstlevel.(task).(['run' run]).plot_params.cond_onset;
%    task_info.durations = options.firstlevel.(task).(['run' run]).plot_params.cond_duration;
%    task_info.precision = 1;
%    trace_info = [];
%
%    roi_fn = roi_anat_fn;
%    trace_info = NFB{sig}(1,:);
%    saveAs_fn = strrep(functional_fn, '_bold.nii', '_anatROItsplot.png');
%    fmrwhy_util_thePlotROIrt(functional_fn, mask_fn, roi_fn, task_info, trace_info, saveAs_fn)
%
%    roi_fn = roi_funcAnat_fn;
%    trace_info = NFB{sig}(2,:);
%    saveAs_fn = strrep(functional_fn, '_bold.nii', '_funcAnatROItsplot.png');
%    fmrwhy_util_thePlotROIrt(functional_fn, mask_fn, roi_fn, task_info, trace_info, saveAs_fn)
%end