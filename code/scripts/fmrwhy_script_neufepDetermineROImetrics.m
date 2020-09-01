% fmrwhy_script_neufepDetermineROIoverlap

% A custom workflow that runs 1st level analysis for all runs of all tasks of specified subjects

% Code steps:
% 1.


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

N_subs = numel(subs);

% -------
% PER TASK and RUN
% -------
% Loop through sessions, tasks, runs, echoes.
for t = 1:numel(tasks)
    task = tasks{t};
    for r = 1:numel(runs)
        run = runs{r};

        disp(['task-' task '_run-' run])

        overlap_summary_fn = fullfile(options.stats_dir, ['sub-all_task-' task '_run-' run '_desc-roiOverlap.tsv']);
        [d, f, e] = fileparts(overlap_summary_fn);
        temp_txt_fn = fullfile(d, [f '.txt']);

        data = nan(N_subs, numel(col_names));

        for s = 1:numel(subs)
            sub = subs{s};
            disp(sub)

            % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
            options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

            % Create text file for saving tmpa values
            tmap_values_fn = fullfile(options.stats_dir, ['sub-' sub '_task-' task '_run-' run '_desc-tmapvalues.tsv']);
            [d1, f1, e1] = fileparts(tmap_values_fn);
            temp_tmapvals_fn = fullfile(d1, [f1 '.txt']);

            % Update workflow params with subject anatomical derivative filenames
            options = fmrwhy_defaults_subAnat(bids_dir, sub, options);
            background_fn = options.rcoregest_anatomical_fn;

            if strcmp(task, 'motor')
                roi_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rleftMotor_roi.nii']);
            else
                roi_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rbilateralAmygdala_roi.nii']);
            end
            nii = nii_tool('load', roi_fn);
            roi_img = fmrwhy_util_createBinaryImg(double(nii.img), 0.1);
            size_anat_roi = sum(roi_img(:));
            data(s, 1) = size_anat_roi;

            % Initialize data for tmap vals per subject
            data_tmapvals = {};
            max_voxels = 0;

            for e = 1:numel(echoes)
                echo = echoes{e};

                % Update directories and filenames
                options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, options.template_echo, options);
                options.sub_dir_stats = fullfile(options.stats_dir, ['sub-' sub]);

                FWE_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-' echo]);
                noFWE_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-' echo '_noFWEp001e20']);

                % Prepare overlap variables and filenames
                consess = options.firstlevel.(task).(['run' run]).contrast_params.consess;
                if numel(consess) > 1
                    k = 3;
                else
                    k = 1;
                end
                str = consess{k}.tcon.name;
                template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
                saveAs_overlap_fn1 = fullfile(FWE_dir_stats, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-' str 'Overlaps' options.roi.(task).desc{1} 'FWE_roi.nii']);
                saveAs_overlap_fn2 = fullfile(noFWE_dir_stats, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-' str 'Overlaps' options.roi.(task).desc{1} 'noFWE_roi.nii']);
                saveAs_montage_fn1 = strrep(saveAs_overlap_fn1, '.nii', '.png');
                saveAs_montage_fn2 = strrep(saveAs_overlap_fn2, '.nii', '.png');

                % Calculate overlap (and create montages)
                fn1 = fullfile(FWE_dir_stats, ['spmT_' sprintf('%04d', k) '_binary_clusters.nii']);
                output = fmrwhy_util_getBinaryOverlap({roi_fn, fn1}, saveAs_overlap_fn1, template_fn, saveAs_montage_fn1, background_fn);
                data(s, 2*e) = output.size_overlap;
                fn2 = fullfile(noFWE_dir_stats, ['spmT_' sprintf('%04d', k) '_binary_clusters.nii']);
                output = fmrwhy_util_getBinaryOverlap({roi_fn, fn2}, saveAs_overlap_fn2, template_fn, saveAs_montage_fn2, background_fn);
                data(s, 2*e+1) = output.size_overlap;

                % Extract T-map values:
                % FWE
                tmap_fn1 = fullfile(FWE_dir_stats, ['spmT_' sprintf('%04d', k) '.nii']);
                cluster_fn1 = fn1;
                nii_tmap = nii_tool('load', tmap_fn1);
                nii_tmap = double(nii_tmap.img(:));
                nii_cluster = nii_tool('load', cluster_fn1);
                binary_cluster = fmrwhy_util_createBinaryImg(double(nii_cluster.img), 0.1);
                I_cluster = find(binary_cluster(:));
                data_tmapvals{2*e-1} = nii_tmap(I_cluster);
                if numel(I_cluster) > max_voxels
                    max_voxels = numel(I_cluster);
                end
                % No FWE
                tmap_fn2 = fullfile(noFWE_dir_stats, ['spmT_' sprintf('%04d', k) '.nii']);
                cluster_fn2 = fn2;
                nii_tmap = nii_tool('load', tmap_fn2);
                nii_tmap = double(nii_tmap.img(:));
                nii_cluster = nii_tool('load', cluster_fn2);
                binary_cluster = fmrwhy_util_createBinaryImg(double(nii_cluster.img), 0.1);
                I_cluster = find(binary_cluster(:));
                data_tmapvals{2*e} = nii_tmap(I_cluster);
                if numel(I_cluster) > max_voxels
                    max_voxels = numel(I_cluster);
                end

            end

            data_tmapvals_mat = nan(max_voxels, numel(data_tmapvals));

            for e = 1:numel(data_tmapvals)
                n = numel(data_tmapvals{e});
                data_tmapvals_mat(1:n, e) = data_tmapvals{e};
            end

            T = array2table(data_tmapvals_mat, 'VariableNames', tmap_colnames);
            writetable(T, temp_tmapvals_fn, 'Delimiter','\t');
            [status, msg, msgID] = movefile(temp_tmapvals_fn, tmap_values_fn);
        end

        data_table = array2table(data,'VariableNames', col_names);
        writetable(data_table, temp_txt_fn, 'Delimiter','\t');
        [status, msg, msgID] = movefile(temp_txt_fn, overlap_summary_fn);
    end
end