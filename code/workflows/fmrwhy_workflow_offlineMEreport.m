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
% STEP 1 -- Load defaults, filenames and parameters
% -------

% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Set subject, sessions
sub = '001';
ses = '';

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

task = 'rest';
run = '1';

% Mask details
masks = fmrwhy_util_loadMasks(bids_dir, sub);
mask_fn = masks.brain_mask_fn;
[p_mask, frm, rg, dim_mask] = fmrwhy_util_readNifti(mask_fn);
mask_img = p_mask.nii.img;
I_mask = find(mask_img(:));

masks_oriented = fmrwhy_util_loadOrientMasks(bids_dir, sub);
mask_fn_oriented = masks_oriented.brain_mask_fn;
[p_mask_oriented, frm_oriented, rg_oriented, dim_mask_oriented] = fmrwhy_util_readNifti(mask_fn_oriented);
mask_img_oriented = p_mask_oriented.nii.img;
I_mask_oriented = find(mask_img_oriented(:));

% -------
% grab template
% -------
template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
options.template_fn = template_fn;

% -------
% Visualise t2star and s0 maps
% -------
t2star_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_t2star.nii']);
[p1, frm1, rg1, dim1] = fmrwhy_util_readOrientNifti(t2star_fn);
t2star_montage = fmrwhy_util_createMontage(p1.nii.img, 9, 1, 'T2star', 'hot', 'off', 'max', [0 200]);
colorbar; % caxis([0 200]);
t2star_png = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_t2star.png']);
if ~exist(t2star_png, 'file')
    print(t2star_montage.f, t2star_png,'-dpng', '-r0')
end

s0_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_s0.nii']);
[p2, frm2, rg2, dim2] = fmrwhy_util_readOrientNifti(s0_fn);
s0_montage = fmrwhy_util_createMontage(p2.nii.img, 9, 1, 'S0', 'parula', 'off', 'max', [0 20000]);
colorbar;
s0_png = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-MEparams_s0.png']);
if ~exist(s0_png, 'file')
    print(s0_montage.f, s0_png,'-dpng', '-r0')
end

% -------
% For all tasks and runs, generate images
% -------
tasks = {'rest', 'motor', 'emotion'};
runs = {'1', '2'};
combined_str = {'Echo 2', 'T2star', 'tSNR', 'TE'};
roi_text = {'', 'left motor cortex', 'bilateral amygdala'};

task_names = {'rest', 'Right finger tapping', 'Hariri task'}

% ROIs
roi_fns = {};
roi_fns{1} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rleftMotor_roi.nii']);
roi_fns{2} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rrightMotor_roi.nii']);
%roi_fns{3} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rleftAmygdala_roi.nii']);
%roi_fns{4} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rrightAmygdala_roi.nii']);
roi_fns{3} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rbilateralAmygdala_roi.nii']);




for t = 1:numel(tasks)
    task = tasks{t};
    for r = 1:numel(runs)
        run = runs{r};
        if t == 1 && r == 1
            continue;
        end
        disp('------------')
        disp('------------')
        disp(['Task: ' task ';  Run: ' run])
        disp('------------')
        disp('------------')

        % Grab filenames for bold, combined, tsnr, smoothed bold/combined
        bold_fns = {};
        bold_fns{1} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-1_desc-rapreproc_bold.nii']);
        bold_fns{2} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_bold.nii']);
        bold_fns{3} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-3_desc-rapreproc_bold.nii']);
        bold_combined_fns{1} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEt2star_bold.nii']);
        bold_combined_fns{2} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEtsnr_bold.nii']);
        bold_combined_fns{3} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEte_bold.nii']);
        tsnr_fns = {};
        tsnr_fns{1} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_tsnr.nii']);
        tsnr_fns{2} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEt2star_tsnr.nii']);
        tsnr_fns{3} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEtsnr_tsnr.nii']);
        tsnr_fns{4} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEte_tsnr.nii']);

        main_fns = [bold_fns{2}, bold_combined_fns];
        smooth_fns = {};
        for i = 1:numel(main_fns)
            smooth_fns{i} = strrep(main_fns{i}, 'desc-', 'desc-s');
        end

        % ---
        % Output 1: Montage of single volume: echo 1, 2, 3, and combined (methods 1, 2, 3)
        % ---
        volume_nr = 5;
        bold_pngs = {};
        for i = 1:numel(bold_fns)
            bold_pngs{i} = strrep(bold_fns{i}, '.nii', '.png');
            if ~exist(bold_pngs{i}, 'file')
                [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(bold_fns{i});
                bold_montage = fmrwhy_util_createMontage(p.nii.img(:,:,:,volume_nr), 9, 1, ['Echo ' num2str(i)], 'gray', 'off', 'max', 0);
                print(bold_montage.f, bold_pngs{i},'-dpng', '-r0');
            end
        end
        combined_pngs = {};
        for i = 1:numel(bold_combined_fns)
            combined_pngs{i} = strrep(bold_combined_fns{i}, '.nii', '.png');
            if ~exist(combined_pngs{i}, 'file')
                [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(bold_combined_fns{i});
                combined_montage = fmrwhy_util_createMontage(p.nii.img(:,:,:,volume_nr), 9, 1, ['Combined - ' combined_str{i+1}], 'gray', 'off', 'max', 0);
                print(combined_montage.f, combined_pngs{i},'-dpng', '-r0');
            end
        end

        % ---
        % Output 2:
        %   - Montage of tSNR: middle echo + combined 1, 2, 3
        %   - Montage of percentage signal change in tsnr: combined 1, 2, 3 vs middle echo
        %   - Raincloud plots of tsnr and percentage signal change in tsnr (tsnr-brain, tsnr-roi, psc,-roi)
        % ---
        tsnr_pngs = {};
        percdiff_pngs = {};
        distr_png = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-tsnrPercdiffRainclouds.png']);
        for i = 1:numel(tsnr_fns)
            tsnr_pngs{i} = strrep(tsnr_fns{i}, '.nii', '.png');
            if i > 1
                percdiff_pngs{i-1} = strrep(tsnr_fns{i}, '_tsnr.nii', '_percdiff.png');
            end
        end

        if strcmp(task, 'rest')
            I_roi_distr = 0;
            roi_img = 0;
        elseif strcmp(task, 'motor')
            overlay_img = {};
            [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(roi_fns{1});
            roi_img = fmrwhy_util_createBinaryImg(p.nii.img, 0.1);
            I_roi_distr = find(roi_img(:));
        else
            roi_img = zeros(dim_mask);
            overlay_img = {};
            [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(roi_fns{3});
            roi_img = fmrwhy_util_createBinaryImg(p.nii.img, 0.1);
%            for i = 3:4
%                [p, frm, rg, dim] = fmrwhy_util_readNifti(roi_fns{i});
%                roi_new = fmrwhy_util_createBinaryImg(p.nii.img, 0.1);
%                roi_img = roi_img | roi_new;
%            end
            I_roi_distr = find(roi_img(:));
        end

        fmrwhy_util_compareTSNR(tsnr_fns, mask_fn, roi_fns, I_roi_distr, roi_text{t} , tsnr_pngs, percdiff_pngs, distr_png);


        % Timeseries plot with task info
        if ~strcmp(task, 'rest')
            str1 = [task_names{t} ' - ' roi_text{t} ' - ' combined_str{1}];
            str2 = [task_names{t} ' - ' roi_text{t} ' - Multi-echo combined ' combined_str{2}];
            str3 = [task_names{t} ' - ' roi_text{t} ' - Multi-echo combined ' combined_str{3}];
            str4 = [task_names{t} ' - ' roi_text{t} ' - Multi-echo combined ' combined_str{4}];
            title_str = {str1, str2, str3, str4};
            for p = 1:numel(smooth_fns)
                functional_fn = smooth_fns{p};
                tsnr_fn = tsnr_fns{p};
            %    str1 = 'Finger tapping - Left motor cortex - Single echo';
            %    str2 = 'Finger tapping - Left motor cortex - Multi-echo combined (method T2star)';
            %    str3 = 'Finger tapping - Left motor cortex - Multi-echo combined (method tSNR)';
            %    str4 = 'Finger tapping - Left motor cortex - Multi-echo combined (method TE)';
                saveAs_fn = strrep(smooth_fns{p}, '_bold.nii', '_tsplot.png');
                task_info.TR = options.firstlevel.(task).sess_params.timing_RT;
                task_info.onsets = options.firstlevel.(task).sess_params.cond_onset;
                task_info.durations = options.firstlevel.(task).sess_params.cond_duration;
                task_info.precision = 1;

                fmrwhy_util_computeROImeasures(functional_fn, roi_img, task_info, tsnr_fn, title_str{p}, saveAs_fn, options)
            end
        end
    end
end

