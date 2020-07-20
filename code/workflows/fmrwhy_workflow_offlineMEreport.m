% A custom workflow that does ...

% Code steps:
% For template task and run (rest, run 1):
%   - T2star montage
%   - S0 montage
% For all other tasks and runs, loop through and generate:
%   - Montage of single volume: echo 1, 2, 3, and combined (methods 1, 2, 3)
%   - Montage of tSNR: middle echo + combined 1, 2, 3
%   - Montage of percentage signal change in tsnr: combined 1, 2, 3 vs middle echo
%   - Raincloud plots of tsnr and percentage signal change in tsnr (tsnr-brain, tsnr-roi, psc,-roi)
%   - Timeseries plot with task info

%--------------------------------------------------------------------------


% -------
% STEP 1: Load defaults, filenames and parameters
% -------

% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);
options.me_dir = fullfile(options.deriv_dir, 'fmrwhy-multiecho');

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Set subject, sessions
sub = '001';
ses = '';

options.sub_dir_me = fullfile(options.me_dir, ['sub-' sub]);
options.func_dir_me = fullfile(options.sub_dir_me, 'func');

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% Plotting settings
rgb_ongray = [255, 115, 236];
rgb_onhot = [148, 239, 255];
rgb_onparula = [255, 115, 236];

% -------
% STEP 2: Grab template data
% -------
task = 'rest';
run = '1';
% Mask details
masks = fmrwhy_util_loadMasks(bids_dir, sub);
mask_fn = masks.brain_mask_fn;
mask_img = masks.brain_mask_3D;
I_mask = masks.brain_mask_I;
masks_oriented = fmrwhy_util_loadOrientMasks(bids_dir, sub);
mask_img_oriented = masks_oriented.brain_mask_3D;
I_mask_oriented = masks_oriented.brain_mask_I;
% Functional volume template
template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
options.template_fn = template_fn;
% ROIs
roi_fns = {};
roi_fns{1} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rleftMotor_roi.nii']);
roi_fns{2} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rbilateralAmygdala_roi.nii']);
compare_roi_txt = {'left motor cortex', 'bilateral amygdala'};
% Grab+load ROI image data; get ROI indices; combine ROI image data into a single overlay image
roi_img = {};
I_roi = {};
overlay_img = zeros(size(mask_img_oriented));
for i = 1:numel(roi_fns)
    [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(roi_fns{i});
    roi_img{i} = fmrwhy_util_createBinaryImg(p.nii.img, 0.1);
    I_roi{i} = find(roi_img{i}(:));
    overlay_img = overlay_img | roi_img{i};
end

% -------
% STEP 3: Visualise t2star and s0 maps
% -------
% t2star
t2star_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_t2star.nii']);
t2star_png = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_t2star.png']);
[p1, frm1, rg1, dim1] = fmrwhy_util_readOrientNifti(t2star_fn);
t2star_img = fmrwhy_util_maskImage(p1.nii.img, mask_img_oriented);
t2star_img(t2star_img>=500) = 0; % TODO, is this fine to do?
t2star_montage = fmrwhy_util_createStatsOverlayMontage(t2star_img, [], overlay_img, 9, 1, '', 'hot', 'off', 'max', [0 200], [], rgb_onhot, true, t2star_png);
%t2star_montage = fmrwhy_util_createMontage(t2star_img, 9, 1, 'T2star', 'hot', 'off', 'max', [0 200]);
%fmrwhy_util_stretchAx(t2star_montage.ax)
%fmrwhy_util_removeTicksAx(t2star_montage.ax)
%colorbar; % caxis([0 200]);




if ~exist(t2star_png, 'file')
    print(t2star_montage.f, t2star_png,'-dpng', '-r0')
end
% S0
s0_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_s0.nii']);
s0_png = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_s0.png']);
[p2, frm2, rg2, dim2] = fmrwhy_util_readOrientNifti(s0_fn);
s0_img = fmrwhy_util_maskImage(p2.nii.img, mask_img_oriented);
s0_montage = fmrwhy_util_createStatsOverlayMontage(s0_img, [], overlay_img, 9, 1, '', 'parula', 'off', 'max', [0 7000], [], rgb_onparula, true, s0_png);
%s0_montage = fmrwhy_util_createMontage(s0_img, 9, 1, 'S0', 'parula', 'off', 'max', [0 20000]);
%fmrwhy_util_stretchAx(s0_montage.ax)
%fmrwhy_util_removeTicksAx(s0_montage.ax)
%colorbar;

if ~exist(s0_png, 'file')
    print(s0_montage.f, s0_png,'-dpng', '-r0')
end

% -------
% STEP 4: For all tasks and runs, generate ME-related images
% -------
tasks = {'rest', 'motor', 'emotion'};
runs = {'1', '2'};
combined_str = {'Echo 2', 'T2star', 'tSNR', 'TE'};
roi_text = {'', 'left motor cortex', 'bilateral amygdala'};
task_names = {'rest', 'Right finger tapping', 'Hariri task'}


for t = 1:numel(tasks)
    task = tasks{t};
    for r = 1:numel(runs)
        run = runs{r};

        % For template task and run
        if strcmp(task, 'rest') == 1 && strcmp(run, '1') == 1
            disp('------------')
            disp(['Task: ' task ';  Run: ' run])
            disp('------------')
            % Grab filenames for bold
            bold_fns = {};
            bold_fns{1} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-1_desc-rapreproc_bold.nii']);
            bold_fns{2} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_bold.nii']);
            bold_fns{3} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-3_desc-rapreproc_bold.nii']);
            % use arbitrary volume number
            volume_nr = 5;
            % Create image outputs for original multi-echo bold data
            bold_pngs = {};
            for i = 1:numel(bold_fns)
                [dir_name, file_name, ext] = fileparts(bold_fns{i});
                bold_pngs{i} = fullfile(options.func_dir_me, [file_name '.png']);
                if ~exist(bold_pngs{i}, 'file')
                    [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(bold_fns{i});
                    bold_img = fmrwhy_util_maskImage(double(p.nii.img(:,:,:,volume_nr)), mask_img_oriented);
                    bold_montage = fmrwhy_util_createStatsOverlayMontage(bold_img, [], overlay_img, 9, 1, '', 'gray', 'off', 'max', [], [], rgb_ongray, false, bold_pngs{i});
                end
            end
            % Create image outputs for original multi-echo tsnr data
            tsnr_fns = {};
            tsnr_pngs = {};
            for i = 1:numel(bold_fns)

                [dir_name, file_name, ext] = fileparts(bold_fns{i});
                tsnr_fns{i} = fullfile(options.func_dir_me, [file_name ext]);
                tsnr_fns{i} = strrep(tsnr_fns{i}, 'bold.nii', 'tsnr.nii');
                tsnr_pngs{i} = strrep(tsnr_fns{i}, '.nii', '.png');
%                [di, fi, ex] = fileparts(tsnr_pngs{i});
%                di
%                fi
%                ex
                if ~exist(tsnr_pngs{i}, 'file')
                    disp(['File exists: ' tsnr_pngs{i}])
                    [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(tsnr_fns{i});
                    tsnr_img = fmrwhy_util_maskImage(double(p.nii.img), mask_img_oriented);
                    tsnr_montage = fmrwhy_util_createStatsOverlayMontage(tsnr_img, [], overlay_img, 9, 1, '', 'hot', 'off', 'max', [0 250], [], rgb_onhot, true, tsnr_pngs{i});
                end
            end

            disp('------------')
            disp(['Skipping Combined files for: Task = ' task ';  Run = ' run])
            disp('------------')
            continue;
        end

        disp('------------')
        disp(['Task: ' task ';  Run: ' run])
        disp('------------')

        % -------
        % STEP 4.1: Single volume outputs
        % -------
        % Output:
        %   - Montages of single volumes of: echo 1, 2, 3
        %   - Montages of single combined volumes of combined timeseries using methods: 1, 2, 3
        % -------
        % Grab filenames for bold and combined
        bold_fns = {};
        bold_fns{1} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-1_desc-rapreproc_bold.nii']);
        bold_fns{2} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_bold.nii']);
        bold_fns{3} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-3_desc-rapreproc_bold.nii']);
        bold_combined_fns{1} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEt2star_bold.nii']);
        bold_combined_fns{2} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEtsnr_bold.nii']);
        bold_combined_fns{3} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEte_bold.nii']);
        main_fns = [bold_fns{2}, bold_combined_fns];
        % use arbitrary volume number
        volume_nr = 5;
        % Create image outputs for original multi-echo data
        bold_pngs = {};
        for i = 1:numel(bold_fns)
            [dir_name, file_name, ext] = fileparts(bold_fns{i});
            bold_pngs{i} = fullfile(options.func_dir_me, [file_name '.png']);
            if ~exist(bold_pngs{i}, 'file')
                [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(bold_fns{i});
                bold_img = fmrwhy_util_maskImage(double(p.nii.img(:,:,:,volume_nr)), mask_img_oriented);
                bold_montage = fmrwhy_util_createStatsOverlayMontage(bold_img, [], overlay_img, 9, 1, '', 'gray', 'off', 'max', [], [], rgb_ongray, false, bold_pngs{i});
%                [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(bold_fns{i});
%                bold_img = fmrwhy_util_maskImage(p.nii.img(:,:,:,volume_nr), mask_img_oriented);
%                bold_montage = fmrwhy_util_createMontage(bold_img, 9, 1, ['Echo ' num2str(i)], 'gray', 'off', 'max', 0);
%                fmrwhy_util_stretchAx(bold_montage.ax)
%                fmrwhy_util_removeTicksAx(bold_montage.ax)
%                print(bold_montage.f, bold_pngs{i},'-dpng', '-r0');
            end
        end
        % Create image outputs for combined multi-echo data
        combined_pngs = {};
        for i = 1:numel(bold_combined_fns)
            combined_pngs{i} = strrep(bold_combined_fns{i}, '.nii', '.png');
            if ~exist(combined_pngs{i}, 'file')
                [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(bold_combined_fns{i});
                combined_img = fmrwhy_util_maskImage(double(p.nii.img(:,:,:,volume_nr)), mask_img_oriented);
                combined_montage = fmrwhy_util_createStatsOverlayMontage(combined_img, [], overlay_img, 9, 1, '', 'gray', 'off', 'max', [], [], rgb_ongray, false, combined_pngs{i});

%                [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(bold_combined_fns{i});
%                combined_img = fmrwhy_util_maskImage(p.nii.img(:,:,:,volume_nr), mask_img_oriented);
%                combined_montage = fmrwhy_util_createMontage(combined_img, 9, 1, ['Combined - ' combined_str{i+1}], 'gray', 'off', 'max', 0);
%                fmrwhy_util_stretchAx(combined_montage.ax)
%                fmrwhy_util_removeTicksAx(combined_montage.ax)
%                print(combined_montage.f, combined_pngs{i},'-dpng', '-r0');
            end
        end

        % -------
        % STEP 4.2: tSNR outputs
        % -------
        % Output:
        %   - Montage of tSNR: middle echo and combined 1, 2, 3
        %   - Montage of percentage signal change in tsnr: combined 1, 2, 3 vs middle echo
        %   - Raincloud plots of tsnr and percentage signal change in tsnr (tsnr-brain, tsnr-roi, psc,-roi)
        % -------
        % Grab filenames for tsnr files
        tsnr_fns = {};
        tsnr_fns{1} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_tsnr.nii']);
        tsnr_fns{2} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEt2star_tsnr.nii']);
        tsnr_fns{3} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEtsnr_tsnr.nii']);
        tsnr_fns{4} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEte_tsnr.nii']);
        tsnr_pngs = {};
        percdiff_pngs = {};
        distr_png = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-tsnrPercdiffRainclouds.png']);
        for i = 1:numel(tsnr_fns)
            tsnr_pngs{i} = strrep(tsnr_fns{i}, '.nii', '.png');
            if i > 1
                percdiff_pngs{i-1} = strrep(tsnr_fns{i}, '_tsnr.nii', '_percdiff.png');
            end
        end

        fmrwhy_util_compareTSNR(tsnr_fns, mask_fn, roi_fns, compare_roi_txt , tsnr_pngs, percdiff_pngs, distr_png);


        % -------
        % STEP 4.3: ROI Timeseries plot outputs
        % -------
        smooth_fns = {};
        for i = 1:numel(main_fns)
            smooth_fns{i} = strrep(main_fns{i}, 'desc-', 'desc-s');
        end
        if ~strcmp(task, 'rest')
            str1 = [task_names{t} ' - ' roi_text{t} ' - ' combined_str{1}];
            str2 = [task_names{t} ' - ' roi_text{t} ' - Multi-echo combined ' combined_str{2}];
            str3 = [task_names{t} ' - ' roi_text{t} ' - Multi-echo combined ' combined_str{3}];
            str4 = [task_names{t} ' - ' roi_text{t} ' - Multi-echo combined ' combined_str{4}];
            title_str = {str1, str2, str3, str4};

            if strcmp(task, 'motor')
                roi_fn = roi_fns{1};
            else
                roi_fn = roi_fns{2};
            end

            for p = 1:numel(smooth_fns)
                functional_fn = smooth_fns{p};
                tsnr_fn = tsnr_fns{p};
                if p == 1
                    [dir_name, file_name, ext] = fileparts(smooth_fns{p});
                    tmp_fn = fullfile(options.func_dir_me, [file_name ext]);
                    saveAs_fn = strrep(tmp_fn, '_bold.nii', '_tsplot.png');
                else
                    saveAs_fn = strrep(smooth_fns{p}, '_bold.nii', '_tsplot.png');
                end
                task_info.TR = options.firstlevel.(task).(['run' run]).sess_params.timing_RT;
                task_info.onsets = options.firstlevel.(task).(['run' run]).plot_params.cond_onset;
                task_info.durations = options.firstlevel.(task).(['run' run]).plot_params.cond_duration;
                task_info.precision = 1;
                if ~exist(saveAs_fn, 'file')
                    trace_info = [];
                    fmrwhy_util_thePlotROI(functional_fn, mask_fn, roi_fn, task_info, trace_info, saveAs_fn)
                else
                    disp(['File already exists: ' saveAs_fn])
                end
            end
        end
    end
end

