% fmrwhy_script_neufepDetermineROIoverlap

% This script performs the following steps for each run of each task (i.e. * 4):
% 1. Per subject, determine the overlap between the functional task clusters from each echo-timeseries and the anatomical task ROI, write all to file: e.g. 'sub-all_task-motor_run-1_desc-roiOverlap.tsv'
% 2.


%--------------------------------------------------------------------------


% -------
% STEP 0.1 -- Load defaults, filenames and parameters
% -------

% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
%bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';
%bids_dir = '/Users/jheunis/Desktop/NEUFEPME_data_BIDS';
bids_dir = '/Volumes/TSM/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Loop through subjects, sessions, tasks, runs, etc
%subs = {'016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
subs = {'001', '002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
%subs = {'001'};
ses = '';
%tasks = {'motor', 'emotion'};
%runs = {'1', '2'};
%echoes = {'2', 'combinedMEtsnr', 'combinedMEt2star', 'combinedMEte'};

tasks = {'motor', 'emotion'};
%tasks = {'motor'};
runs = {'1', '2'};
%runs = {'1'};
echoes = {'2', 'combinedMEtsnr', 'combinedMEt2star', 'combinedMEte', 'combinedMEt2starFIT', 't2starFIT'};
%echoes = {'2'};
rgb_onhot = [148, 239, 255];

col_names = {'anat_roi', 'echo2_FWE', 'echo2_noFWE', 'combTSNR_FWE', 'combTSNR_noFWE', 'combT2STAR_FWE', 'combT2STAR_noFWE', 'combTE_FWE', 'combTE_noFWE', 'combT2STARfit_FWE', 'combT2STARfit_noFWE', 'T2STARfit_FWE', 'T2STARfit_noFWE'};

tmap_colnames = {'echo2_FWE', 'echo2_noFWE', 'combTSNR_FWE', 'combTSNR_noFWE', 'combT2STAR_FWE', 'combT2STAR_noFWE', 'combTE_FWE', 'combTE_noFWE', 'combT2STARfit_FWE', 'combT2STARfit_noFWE', 'T2STARfit_FWE', 'T2STARfit_noFWE'};

echo_colnames = {'echo2', 'combTSNR', 'combT2STAR', 'combTE', 'combT2STARfit', 'T2STARfit'};
cluster_colnames = {'FWE', 'noFWE', 'anatROI', 'fweAND', 'fweOR'};

N_subs = numel(subs);

% -------
% PER TASK and RUN
% -------
% Loop through tasks, runs
for t = 1:numel(tasks)
    task = tasks{t};
    for r = 1:numel(runs)
        run = runs{r};

        disp(['task-' task '_run-' run])

        % Initialize file and matrix for roi overlap data
        overlap_summary_fn = fullfile(options.stats_dir, ['sub-all_task-' task '_run-' run '_desc-roiOverlap.tsv']);
        [d, f, e] = fileparts(overlap_summary_fn);
        temp_txt_fn = fullfile(d, [f '.txt']);
        roi_overlap_data = nan(N_subs, numel(col_names));
        % Initialize files and matrices for mean and peak data
        meanCvalue_mat = [];
        meanTvalue_mat = [];
        meanPSCvalue_mat = [];
        peakCvalue_mat = [];
        peakTvalue_mat = [];
        peakPSCvalue_mat = [];


        for s = 1:numel(subs)
            sub = subs{s};
            disp(sub)

            % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
            options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

            % Create text files for saving tmap, effect size (con maps), and PSC (map+time series) values
            tmap_values_fn = fullfile(options.stats_dir, ['sub-' sub '_task-' task '_run-' run '_desc-tmapvalues.tsv']);
            cmap_values_fn = fullfile(options.stats_dir, ['sub-' sub '_task-' task '_run-' run '_desc-cmapvalues.tsv']);
            psc_values_fn = fullfile(options.stats_dir, ['sub-' sub '_task-' task '_run-' run '_desc-PSCvalues.tsv']);
            psc_timeseries_fn = fullfile(options.stats_dir, ['sub-' sub '_task-' task '_run-' run '_desc-PSCtimeseries.tsv']);

            echo_cluster_colnames = {};

            % Update workflow params with subject anatomical derivative filenames
            options = fmrwhy_defaults_subAnat(bids_dir, sub, options);
            background_fn = options.rcoregest_anatomical_fn;
            [pbackground, ~, ~, ~] = fmrwhy_util_readOrientNifti(background_fn);
            background_img = double(pbackground.nii.img);

            % Grab ROI fn
            if strcmp(task, 'motor')
                roi_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rleftMotor_roi.nii']);
            else
                roi_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rbilateralAmygdala_roi.nii']);
            end

            % Load task ROI and set first column of data (for the subject) as the sum of voxels in ROI
            % TODO: there is a difference in voxel values when loading the resliced ROI files using nii_tool versus fmrwhy_util_readOrientNifti
            % TODO: this issue might result from the roi files not being binary, but rather some probabalistic value between 0 and 1. When thresholding these images, and then summing the 1s, this yeilds different outcomes! Need to fix!
            % TODO: the nii_tool result seems to be the correct one, since this yields a number of voxels in the ROI that is plausible, the other option yields a number that seems too low.
            nii = nii_tool('load', roi_fn);
            roi_img = fmrwhy_util_createBinaryImg(double(nii.img), 0.1);
            size_anat_roi = sum(roi_img(:));
            roi_overlap_data(s, 1) = size_anat_roi;

            % Grab mask fn
            mask_fn = options.brain_mask_fn;

            % Grab template fn
            template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);

            % Initialize data for tmap, effect size (contrast), and PSC (map+time series) values, per subject
            data_tmapvals = {};
            data_cmapvals = {};
            data_pscvals = {};
            data_pscts = nan(options.Nscans, numel(echoes)*numel(cluster_colnames));

            % Initialize filenames for BOLD time series to extract PSC time series
            options.func_dir_me = fullfile(options.me_dir, ['sub-' sub], 'func' );
            functional_fns = {};
            functional_fns{1} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-srapreproc_bold.nii']);
            functional_fns{2} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-scombinedMEt2star_bold.nii']);
            functional_fns{3} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-scombinedMEtsnr_bold.nii']);
            functional_fns{4} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-scombinedMEte_bold.nii']);
            functional_fns{5} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-scombinedMEt2starFIT_bold.nii']);
            functional_fns{6} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-st2starFIT_bold.nii']);

            % Variable to store max voxels, to be used to size the resulting data matrix appropriately
            max_voxels = 0;

            % echo_cluster_colnames for eventual text file writing
            echo_cluster_colnames = {};

            %----------
            %----------
            % Loop through echo options and create task activation maps overlap images (logical OR, logical AND)
            %----------
            % Based on FWE thresholds
            %----------
            binary_fns = {};
            binary_imgs = {};

            binary_OR_fn = fullfile(options.sub_dir_stats, ['sub-' sub '_task-' task '_run-' run '_desc-allEchoTaskClustersOR.nii']);
            binary_AND_fn = fullfile(options.sub_dir_stats, ['sub-' sub '_task-' task '_run-' run '_desc-allEchoTaskClustersAND.nii']);
            binary_SUM_fn = fullfile(options.sub_dir_stats, ['sub-' sub '_task-' task '_run-' run '_desc-allEchoTaskClustersSUM.nii']);
            binary_OR_fn_noFWE = fullfile(options.sub_dir_stats, ['sub-' sub '_task-' task '_run-' run '_desc-allEchoTaskClustersOR_noFWE.nii']);
            binary_AND_fn_noFWE = fullfile(options.sub_dir_stats, ['sub-' sub '_task-' task '_run-' run '_desc-allEchoTaskClustersAND_noFWE.nii']);
            binary_SUM_fn_noFWE = fullfile(options.sub_dir_stats, ['sub-' sub '_task-' task '_run-' run '_desc-allEchoTaskClustersSUM_noFWE.nii']);

            for e = 1:numel(echoes)
                echo = echoes{e};

                % Update directories and filenames
                options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, options.template_echo, options);
                options.sub_dir_stats = fullfile(options.stats_dir, ['sub-' sub]);
                FWE_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-' echo]);
                noFWE_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-' echo '_noFWEp001e20']);

                % Grab correct contrast number
                consess = options.firstlevel.(task).(['run' run]).contrast_params.consess;
                if numel(consess) > 1
                    k = 3;
                else
                    k = 1;
                end
                str = consess{k}.tcon.name;

                % Get image data
                binary_fns{e} = fullfile(FWE_dir_stats, ['spmT_' sprintf('%04d', k) '_binary_clusters.nii']);
                binary_fn_noFWE = fullfile(noFWE_dir_stats, ['spmT_' sprintf('%04d', k) '_binary_clusters.nii']);

                nii = nii_tool('load', binary_fns{e});
                nii_noFWE = nii_tool('load', binary_fn_noFWE);
                % Make sure the image is in binary format
                binary_imgs{e} = fmrwhy_util_createBinaryImg(double(nii.img), 0.1);
                binary_img_noFWE = fmrwhy_util_createBinaryImg(double(nii_noFWE.img), 0.1);

                % Get logical AND and OR of new image with previous images
                if e == 1
                    binary_OR_img = binary_imgs{e};
                    binary_AND_img = binary_imgs{e};
                    sum_img = binary_imgs{e};
                    binary_OR_img_noFWE = binary_img_noFWE;
                    binary_AND_img_noFWE = binary_img_noFWE;
                    sum_img_noFWE = binary_img_noFWE;
                else
                    binary_OR_img = binary_OR_img | binary_imgs{e};
                    binary_AND_img = binary_AND_img & binary_imgs{e};
                    binary_OR_img_noFWE = binary_OR_img_noFWE | binary_img_noFWE;
                    binary_AND_img_noFWE = binary_AND_img_noFWE & binary_img_noFWE;
                end

                sum_img = sum_img + binary_imgs{e};
                sum_img_noFWE = sum_img_noFWE + binary_img_noFWE;
            end

            % Calculate the overlap of cluster images as percentage
            sum_img = sum_img/numel(echoes)*100;
            sum_img_noFWE = sum_img_noFWE/numel(echoes)*100;

            % Save images
            if ~exist(binary_SUM_fn, 'file')
                no_scaling = 1;
                fmrwhy_util_saveNifti(binary_SUM_fn, double(sum_img), template_fn, no_scaling);
                [psum, ~, ~, ~] = fmrwhy_util_readOrientNifti(binary_SUM_fn);
                sum_img = double(psum.nii.img);
                binary_SUM_png = strrep(binary_SUM_fn, '.nii', '.png');
                sum_montage = fmrwhy_util_createStatsOverlayMontage(background_img, sum_img, [], 9, 1, '', 'gray', 'off', 'max', [0 100], 'spring', [], true, binary_SUM_png);

                fmrwhy_util_saveNifti(binary_SUM_fn_noFWE, double(sum_img_noFWE), template_fn, no_scaling);
                [psum_noFWE, ~, ~, ~] = fmrwhy_util_readOrientNifti(binary_SUM_fn_noFWE);
                sum_img_noFWE = double(psum_noFWE.nii.img);
                binary_SUM_png_noFWE = strrep(binary_SUM_fn_noFWE, '.nii', '.png');
                sum_montage_noFWE = fmrwhy_util_createStatsOverlayMontage(background_img, sum_img_noFWE, [], 9, 1, '', 'gray', 'off', 'max', [0 100], 'spring', [], true, binary_SUM_png_noFWE);

            end
            if ~exist(binary_OR_fn, 'file') || ~exist(binary_AND_fn, 'file')
                no_scaling = 1;
                roi_rgbcolors = [148, 239, 255];

                fmrwhy_util_saveNifti(binary_OR_fn, double(binary_OR_img), template_fn, no_scaling);
                fmrwhy_util_saveNifti(binary_AND_fn, double(binary_AND_img), template_fn, no_scaling);

                [por, ~, ~, ~] = fmrwhy_util_readOrientNifti(binary_OR_fn);
                or_img = double(por.nii.img);
                [pand, ~, ~, ~] = fmrwhy_util_readOrientNifti(binary_AND_fn);
                and_img = double(pand.nii.img);
                binary_OR_png = strrep(binary_OR_fn, '.nii', '.png');
                binary_AND_png = strrep(binary_AND_fn, '.nii', '.png');

                overlapmontage = fmrwhy_util_createStatsOverlayMontage(background_img, [], or_img, 9, 1, '', 'gray', 'off', 'max', [], [], roi_rgbcolors, false, binary_OR_png);
                overlapmontage = fmrwhy_util_createStatsOverlayMontage(background_img, [], and_img, 9, 1, '', 'gray', 'off', 'max', [], [], roi_rgbcolors, false, binary_AND_png);

                fmrwhy_util_saveNifti(binary_OR_fn_noFWE, double(binary_OR_img_noFWE), template_fn, no_scaling);
                fmrwhy_util_saveNifti(binary_AND_fn_noFWE, double(binary_AND_img_noFWE), template_fn, no_scaling);

                [por_noFWE, ~, ~, ~] = fmrwhy_util_readOrientNifti(binary_OR_fn_noFWE);
                or_img_noFWE = double(por_noFWE.nii.img);
                [pand_noFWE, ~, ~, ~] = fmrwhy_util_readOrientNifti(binary_AND_fn_noFWE);
                and_img_noFWE = double(pand_noFWE.nii.img);
                binary_OR_png_noFWE = strrep(binary_OR_fn_noFWE, '.nii', '.png');
                binary_AND_png_noFWE = strrep(binary_AND_fn_noFWE, '.nii', '.png');

                overlapmontage = fmrwhy_util_createStatsOverlayMontage(background_img, [], or_img_noFWE, 9, 1, '', 'gray', 'off', 'max', [], [], roi_rgbcolors, false, binary_OR_png_noFWE);
                overlapmontage = fmrwhy_util_createStatsOverlayMontage(background_img, [], and_img_noFWE, 9, 1, '', 'gray', 'off', 'max', [], [], roi_rgbcolors, false, binary_AND_png_noFWE);
            end
            %----------


            %----------
            %----------
            % Loop through echo options and extract data for effect size, PSC, etc, from clusters
            %----------
            for e = 1:numel(echoes)
                echo = echoes{e};

                % Update directories and filenames
                options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, options.template_echo, options);
                options.sub_dir_stats = fullfile(options.stats_dir, ['sub-' sub]);

                FWE_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-' echo]);
                noFWE_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-' echo '_noFWEp001e20']);

                % Grab correct contrast number
                consess = options.firstlevel.(task).(['run' run]).contrast_params.consess;
                if numel(consess) > 1
                    k = 3;
                else
                    k = 1;
                end
                str = consess{k}.tcon.name;

                % Overlap filenames
                saveAs_overlap_fn1 = fullfile(FWE_dir_stats, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-' str 'Overlaps' options.roi.(task).desc{1} 'FWE_roi.nii']);
                saveAs_overlap_fn2 = fullfile(noFWE_dir_stats, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-' str 'Overlaps' options.roi.(task).desc{1} 'noFWE_roi.nii']);
                saveAs_montage_fn1 = strrep(saveAs_overlap_fn1, '.nii', '.png');
                saveAs_montage_fn2 = strrep(saveAs_overlap_fn2, '.nii', '.png');

                % Calculate func+anat overlap (and create montages); only if these files do not yet exist
                cluster_fn1 = fullfile(FWE_dir_stats, ['spmT_' sprintf('%04d', k) '_binary_clusters.nii']);
                if ~exist(saveAs_overlap_fn1, 'file')
                    output = fmrwhy_util_getBinaryOverlap({roi_fn, cluster_fn1}, saveAs_overlap_fn1, template_fn, saveAs_montage_fn1, background_fn);
                else
                    output = fmrwhy_util_getBinaryOverlap({roi_fn, cluster_fn1});
                end
                roi_overlap_data(s, 2*e) = output.size_overlap;
                cluster_fn2 = fullfile(noFWE_dir_stats, ['spmT_' sprintf('%04d', k) '_binary_clusters.nii']);
                if ~exist(saveAs_overlap_fn1, 'file')
                    output = fmrwhy_util_getBinaryOverlap({roi_fn, cluster_fn2}, saveAs_overlap_fn2, template_fn, saveAs_montage_fn2, background_fn);
                else
                    output = fmrwhy_util_getBinaryOverlap({roi_fn, cluster_fn2});
                end
                roi_overlap_data(s, 2*e+1) = output.size_overlap;

                %----------
                %----------
                % Generate PSC maps: TODO: make this into a function
                %----------
                % Only done in FWE directory
                %----------
                % Calculate scaling factor for PSC: see (1) http://www.sbirc.ed.ac.uk/cyril/bold_percentage/BOLD_percentage.html and (2) https://www.frontiersin.org/articles/10.3389/fnins.2014.00001/full
                spm_fn = fullfile(FWE_dir_stats, 'SPM.mat');
                spm = load(spm_fn);
                block_length = 20; % seconds
                ntime = block_length*(1/spm.SPM.xBF.dt);
                reference_block =  conv(ones(1,ntime),spm.SPM.xBF.bf(:,1))';
                scale_factor = max(reference_block);

                % Grab beta map of constant regressor
                beta_constant_fn = fullfile(FWE_dir_stats, ['beta_' sprintf('%04d', 19) '.nii']);
                if k == 3
                    beta_constant_fn = fullfile(FWE_dir_stats, ['beta_' sprintf('%04d', 20) '.nii']);
                end
                nii = nii_tool('load', beta_constant_fn);
                constant_vals = double(nii.img);

                % Load con maps, calculate PSC, create nifti
                PSC = {};
                CON = {};
                PSC_2d = {};
                for i = 1:k
                    con_fn{i} = fullfile(FWE_dir_stats, ['con_' sprintf('%04d', i) '.nii']);
                    nii = nii_tool('load', con_fn{i});
                    CON{i} = double(nii.img);
                    PSC{i} = CON{i} * SF ./ constant_vals * 100;
                    psc_img_fn = fullfile(FWE_dir_stats, ['PSC_' sprintf('%04d', i) '.nii']);
                    no_scaling = 1;
                    fmrwhy_util_saveNifti(psc_img_fn, PSC{i}, con_fn{i}, no_scaling)
                    PSC_2d{i} = PSC{i}(:);
                end
                %----------
                %----------


                %----------
                %----------
                % Define clusters
                %----------
                cluster_fns = {cluster_fn1, cluster_fn2, roi_fn, binary_AND_fn, binary_OR_fn};
                % Load cluster image 1 - task activation, FWE
                nii_cluster1 = nii_tool('load', cluster_fn1);
                binary_cluster1 = fmrwhy_util_createBinaryImg(double(nii_cluster1.img), 0.1);
                % Load cluster image 2 - task activation, no FWE
                nii_cluster2 = nii_tool('load', cluster_fn2);
                binary_cluster2 = fmrwhy_util_createBinaryImg(double(nii_cluster2.img), 0.1);
                % Cluster images 3,4,5 already loaded
                cluster_imgs = {binary_cluster1, binary_cluster2, roi_img, binary_AND_img, binary_OR_img}; % the AND/OR images are from FWE thresholded stat maps; would look different from noFWE maps
                %----------

                %----------
                %----------
                % Load data: effect size, t-map, psc values, etc
                %----------
                % Effect sizes (contrast maps)
                cmap_fn1 = fullfile(FWE_dir_stats, ['con_' sprintf('%04d', k) '.nii']);
                nii_cmap = nii_tool('load', cmap_fn1);
                img_cmap = double(nii_cmap.img(:));
                % T-maps
                tmap_fn1 = fullfile(FWE_dir_stats, ['spmT_' sprintf('%04d', k) '.nii']);
                nii_tmap = nii_tool('load', tmap_fn1);
                img_tmap = double(nii_tmap.img(:));

                %----------
                %----------
                % Loop through clusters and extract effect size, t-map, psc values from various clusters:
                % TODO: should tsnr also be extracted here?
                %----------
                Nc = numel(cluster_colnames);
                for c = 1:Nc

                    % add colnames
                    echo_cluster_colnames{Nc*e-Nc+c} = [echo_colnames{e} '_' cluster_colnames{c}];
                    % set max_voxels
                    I_cluster = find(cluster_imgs{c}(:));
                    if numel(I_cluster) > max_voxels
                        max_voxels = numel(I_cluster);
                    end
                    % extract data
                    data_cmapvals{Nc*e-Nc+c} = img_cmap(I_cluster);
                    data_tmapvals{Nc*e-Nc+c} = img_tmap(I_cluster);
                    data_pscvals{Nc*e-Nc+c} = PSC_2d{k}(I_cluster);
                    data_pscts(:, Nc*e-Nc+c) = fmrwhy_util_getROItimeseries(functional_fns{e}, mask_fn, cluster_fns{c});
                end
            end

            % Create matrices with correct sizes based on max_voxels
            data_cmapvals_mat = nan(max_voxels, numel(data_cmapvals));
            data_tmapvals_mat = nan(max_voxels, numel(data_tmapvals));
            data_pscvals_mat = nan(max_voxels, numel(data_pscvals));

            % write data from cell arays to matrices
            for j = 1:numel(data_cmapvals)
                n = numel(data_cmapvals{j});
                data_cmapvals_mat(1:n, j) = data_cmapvals{j};
                data_tmapvals_mat(1:n, j) = data_tmapvals{j};
                data_pscvals_mat(1:n, j) = data_pscvals{j};
            end

            % Create temp text files
            [d1, f1, e1] = fileparts(cmap_values_fn);
            [d2, f2, e2] = fileparts(tmap_values_fn);
            [d3, f3, e3] = fileparts(psc_values_fn);
            [d4, f4, e4] = fileparts(psc_timeseries_fn);
            temp_cmapvals_fn = fullfile(d1, [f1 '.txt']);
            temp_tmapvals_fn = fullfile(d2, [f2 '.txt']);
            temp_pscvals_fn = fullfile(d3, [f3 '.txt']);
            temp_pscts_fn = fullfile(d4, [f4 '.txt']);

            % Write effect sizes to tsv file
            C = array2table(data_cmapvals_mat, 'VariableNames', echo_cluster_colnames);
            writetable(C, temp_cmapvals_fn, 'Delimiter','\t');
            [status, msg, msgID] = movefile(temp_cmapvals_fn, cmap_values_fn);

            % Write T-values to tsv file
            T = array2table(data_tmapvals_mat, 'VariableNames', echo_cluster_colnames);
            writetable(T, temp_tmapvals_fn, 'Delimiter','\t');
            [status, msg, msgID] = movefile(temp_tmapvals_fn, tmap_values_fn);

            % Write PSC values to tsv file
            P = array2table(data_pscvals_mat, 'VariableNames', echo_cluster_colnames);
            writetable(T, temp_pscvals_fn, 'Delimiter','\t');
            [status, msg, msgID] = movefile(temp_pscvals_fn, psc_values_fn);

            % Write PSC timeseries to tsv file
            Pts = array2table(data_pscts, 'VariableNames', echo_cluster_colnames);
            writetable(Pts, temp_pscts_fn, 'Delimiter','\t');
            [status, msg, msgID] = movefile(temp_pscts_fn, psc_timeseries_fn);

            % Calculate mean and peak values and add to matrices
            meanCvalue_mat = [meanCvalue_mat; nanmean(data_cmapvals_mat, 1);];
            meanTvalue_mat = [meanTvalue_mat; nanmean(data_tmapvals_mat, 1);];
            meanPSCvalue_mat = [meanPSCvalue_mat; nanmean(data_pscvals_mat, 1);];
            peakCvalue_mat = [peakCvalue_mat; max(data_cmapvals_mat,[],1);];
            peakTvalue_mat = [peakTvalue_mat; max(data_tmapvals_mat,[],1);];
            peakPSCvalue_mat = [peakPSCvalue_mat; max(data_pscvals_mat,[],1);];

        end

        % Write ROI overlap data to tsv file
        data_table = array2table(data,'VariableNames', col_names);
        writetable(data_table, temp_txt_fn, 'Delimiter','\t');
        [status, msg, msgID] = movefile(temp_txt_fn, overlap_summary_fn);

        % Write mean/peak data to tsv file
        meanCvalue_fn = fullfile(options.stats_dir, ['sub-all_task-' task '_run-' run '_desc-meanCvalues.tsv']);
        meanTvalue_fn = fullfile(options.stats_dir, ['sub-all_task-' task '_run-' run '_desc-meanTvalues.tsv']);
        meanPSCvalue_fn = fullfile(options.stats_dir, ['sub-all_task-' task '_run-' run '_desc-meanPSCCvalues.tsv']);
        peakCvalue_fn = fullfile(options.stats_dir, ['sub-all_task-' task '_run-' run '_desc-peakCvalues.tsv']);
        peakTvalue_fn = fullfile(options.stats_dir, ['sub-all_task-' task '_run-' run '_desc-peakTvalues.tsv']);
        peakPSCvalue_fn = fullfile(options.stats_dir, ['sub-all_task-' task '_run-' run '_desc-peakPSCvalues.tsv']);

        meanPeak_mats = {meanCvalue_mat, meanTvalue_mat, meanPSCvalue_mat, peakCvalue_mat, peakTvalue_mat, peakPSCvalue_mat};
        meanPeak_fns = {meanCvalue_fn, meanTvalue_fn, meanPSCvalue_fn, peakCvalue_fn, peakTvalue_fn, peakPSCvalue_fn};

        for v = 1:numel(meanPeak_mats)
            [d, f, e] = fileparts(meanPeak_fns{v});
            temp_txt_fn = fullfile(d, [f '.txt']);
            data_table = array2table(meanPeak_mats{v},'VariableNames', echo_cluster_colnames);
            writetable(data_table, temp_txt_fn, 'Delimiter','\t');
            [status, msg, msgID] = movefile(temp_txt_fn, meanPeak_fns{v});
        end
    end
end