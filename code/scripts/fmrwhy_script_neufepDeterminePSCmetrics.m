% fmrwhy_script_neufepDeterminePSCmetrics

% A custom workflow that runs 1st level analysis for all runs of all tasks of specified subjects

% Code steps:
% 1.


%--------------------------------------------------------------------------

% TODO: build this pipeline into the fmrwhy_script_neufepDetermineROImetrics script
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
tasks = {'motor', 'emotion'};
%tasks = {'motor'};
runs = {'1', '2'};
%runs = {'1'};
echoes = {'2', 'combinedMEtsnr', 'combinedMEt2star', 'combinedMEte', 'combinedMEt2starFIT', 't2starFIT'};
%echoes = {'2'};
rgb_onhot = [148, 239, 255];

col_names = {'anat_roi', 'echo2_FWE', 'echo2_noFWE', 'combTSNR_FWE', 'combTSNR_noFWE', 'combT2STAR_FWE', 'combT2STAR_noFWE', 'combTE_FWE', 'combTE_noFWE', 'combT2STARfit_FWE', 'combT2STARfit_noFWE', 'T2STARfit_FWE', 'T2STARfit_noFWE'};

tmap_colnames = {'echo2_FWE', 'echo2_noFWE', 'combTSNR_FWE', 'combTSNR_noFWE', 'combT2STAR_FWE', 'combT2STAR_noFWE', 'combTE_FWE', 'combTE_noFWE', 'combT2STARfit_FWE', 'combT2STARfit_noFWE', 'T2STARfit_FWE', 'T2STARfit_noFWE'};
psc_ts_colnames = {'echo2', 'combTSNR', 'combT2STAR', 'combTE', 'combT2STARfit', 'T2STARfit'};

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

        for s = 1:numel(subs)
            sub = subs{s};
            disp(sub)

            % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
            options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

            % Create text file for saving PSC values
            psc_values_fn = fullfile(options.stats_dir, ['sub-' sub '_task-' task '_run-' run '_desc-PSCvalues.tsv']);
            psc_timeseries_fn = fullfile(options.stats_dir, ['sub-' sub '_task-' task '_run-' run '_desc-PSCtimeseries.tsv']);

            [d1, f1, e1] = fileparts(psc_values_fn);
            [d2, f2, e2] = fileparts(psc_timeseries_fn);
            temp_pscvals_fn = fullfile(d1, [f1 '.txt']);
            temp_pscts_fn = fullfile(d2, [f2 '.txt']);

            % Update workflow params with subject anatomical derivative filenames
            options = fmrwhy_defaults_subAnat(bids_dir, sub, options);
            background_fn = options.rcoregest_anatomical_fn;

            % Grab mask fn
            mask_fn = options.brain_mask_fn;

            % Grab ROIs
            if strcmp(task, 'motor')
                roi_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rleftMotor_roi.nii']);
            else
                roi_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rbilateralAmygdala_roi.nii']);
            end
%            nii = nii_tool('load', roi_fn);
%            roi_img = fmrwhy_util_createBinaryImg(double(nii.img), 0.1);
%            size_anat_roi = sum(roi_img(:));

            % Initialize data per subject
            data_pscvals = {};
            data_pscts = nan(options.Nscans, numel(echoes));
            max_voxels = 0;
            options.func_dir_me = fullfile(options.me_dir, ['sub-' sub], 'func' );
            functional_fns = {};
            functional_fns{1} = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-srapreproc_bold.nii']);
            functional_fns{2} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-scombinedMEt2star_bold.nii']);
            functional_fns{3} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-scombinedMEtsnr_bold.nii']);
            functional_fns{4} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-scombinedMEte_bold.nii']);
            functional_fns{5} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-scombinedMEt2starFIT_bold.nii']);
            functional_fns{6} = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-st2starFIT_bold.nii']);

            for e = 1:numel(echoes)
                echo = echoes{e};

                % Update directories and filenames
                options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, options.template_echo, options);
                options.sub_dir_stats = fullfile(options.stats_dir, ['sub-' sub]);

                FWE_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-' echo]);
                noFWE_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-' echo '_noFWEp001e20']);

                % Prepare variables and filenames
                consess = options.firstlevel.(task).(['run' run]).contrast_params.consess;
                if numel(consess) > 1
                    k = 3;
                else
                    k = 1;
                end

%                str = consess{k}.tcon.name;
%                template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);

                % Get binary clusters, and get max
                cluster_fn1 = fullfile(FWE_dir_stats, ['spmT_' sprintf('%04d', k) '_binary_clusters.nii']);
                cluster_fn2 = fullfile(noFWE_dir_stats, ['spmT_' sprintf('%04d', k) '_binary_clusters.nii']);
                nii_cluster1 = nii_tool('load', cluster_fn1);
                binary_cluster1 = fmrwhy_util_createBinaryImg(double(nii_cluster1.img), 0.1);
                I_cluster1 = find(binary_cluster1(:));
                if numel(I_cluster1) > max_voxels
                    max_voxels = numel(I_cluster1);
                end
                nii_cluster2 = nii_tool('load', cluster_fn2);
                binary_cluster2 = fmrwhy_util_createBinaryImg(double(nii_cluster2.img), 0.1);
                I_cluster2 = find(binary_cluster2(:));
                if numel(I_cluster2) > max_voxels
                    max_voxels = numel(I_cluster2);
                end

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

                % Load con maps, calculate PSC, create nifti, save TSV
                PSC = {};
                CON = {};
                PSC_2d = {};
                for i = 1:k
                    con_fn{i} = fullfile(FWE_dir_stats, ['con_' sprintf('%04d', i) '.nii']);
                    nii = nii_tool('load', con_fn{i});
                    CON{i} = double(nii.img);
                    PSC{i} = CON{i} * scale_factor ./ constant_vals * 100;
                    psc_img_fn = fullfile(FWE_dir_stats, ['PSC_' sprintf('%04d', i) '.nii']);
                    no_scaling = 1;
                    fmrwhy_util_saveNifti(psc_img_fn, PSC{i}, con_fn{i}, no_scaling)
                    PSC_2d{i} = PSC{i}(:);
                end

                % FWE
                data_pscvals{2*e-1} = PSC_2d{k}(I_cluster1);
                % noFWE
                data_pscvals{2*e} = PSC_2d{k}(I_cluster2);

                % ROI time series
                data_pscts(:, e) = fmrwhy_util_getROItimeseries(functional_fns{e}, mask_fn, roi_fn);
            end

            % Write PSC values to tsv
            data_pscvals_mat = nan(max_voxels, numel(data_pscvals));
            for e = 1:numel(data_pscvals)
                n = numel(data_pscvals{e});
                data_pscvals_mat(1:n, e) = data_pscvals{e};
            end
            T2 = array2table(data_pscvals_mat, 'VariableNames', tmap_colnames);
            writetable(T2, temp_pscvals_fn, 'Delimiter','\t');
            [status, msg, msgID] = movefile(temp_pscvals_fn, psc_values_fn);

            % Write PSC timeseries to tsv
            T3 = array2table(data_pscts, 'VariableNames', psc_ts_colnames);
            writetable(T3, temp_pscts_fn, 'Delimiter','\t');
            [status, msg, msgID] = movefile(temp_pscts_fn, psc_timeseries_fn);
        end
    end
end