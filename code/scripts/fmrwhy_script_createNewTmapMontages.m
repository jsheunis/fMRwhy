% fmrwhy_script_createNewTmapMontages

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
bids_dir = '/Users/jheunis/Desktop/NEUFEPME_data_BIDS';

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
echoes = {'2', 'combinedMEtsnr', 'combinedMEt2star', 'combinedMEte'};
%echoes = {'2'};
rgb_onhot = [148, 239, 255];

tic;
for s = 1:numel(subs)
    sub = subs{s};
    % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

    roi_fns = {};
    roi_fns{1} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rleftMotor_roi.nii']);
    roi_fns{2} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rbilateralAmygdala_roi.nii']);
    % Grab+load ROI image data; get ROI indices; combine ROI image data into a single overlay image
    roi_img = {};
    I_roi = {};
    for i = 1:numel(roi_fns)
        [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(roi_fns{i});
        roi_img{i} = fmrwhy_util_createBinaryImg(p.nii.img, 0.1);
        I_roi{i} = find(roi_img{i}(:));

        if i == 1
            combined_roi_img = zeros(size(roi_img{i}));
        end

        combined_roi_img = combined_roi_img | roi_img{i};
    end

    % -------
    % PER TASK and RUN
    % -------
    % Loop through sessions, tasks, runs, echoes.

    toTransform_fns = {};
    saveAs_transform_fns = {};
    for t = 1:numel(tasks)
        task = tasks{t};
        for r = 1:numel(runs)
            run = runs{r};
            for e = 1:numel(echoes)
                echo = echoes{e};

                txt = ['Running analysis for: sub-' sub '_task-' task '_run-' run '_echo-' echo];
                disp('------------------')
                disp(txt)
                disp('------------------')
                options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, options.template_echo, options);
                options.sub_dir_stats = fullfile(options.stats_dir, ['sub-' sub]);

                run_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-' echo]);
                run_dir_stats_new = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-' echo '_noFWEp001e20']);

                consess = options.firstlevel.(task).(['run' run]).contrast_params.consess;
                background_fn = options.rcoregest_anatomical_fn;
                [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(background_fn);
                background_img = p.nii.img;

                for k = 1:numel(consess)
                    tmap_fn = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', k) '.nii']);
                    tmap_clusters_fn = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', k) '_binary_clusters.nii']);
                    [ptmap, ~, ~, ~] = fmrwhy_util_readOrientNifti(tmap_fn);
                    [ptmapc, ~, ~, ~] = fmrwhy_util_readOrientNifti(tmap_clusters_fn);
                    stats_img = fmrwhy_util_maskImage(double(ptmap.nii.img), double(ptmapc.nii.img));
                    str = consess{k}.tcon.name;
                    saveAs_fn = fullfile(run_dir_stats_new, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-' str '_threshtmap.png']);
                    overlaymontage = fmrwhy_util_createStatsOverlayMontage(p.nii.img, stats_img, combined_roi_img, 9, 1, '', 'gray', 'off', 'max', [], 'hot', rgb_onhot, true, saveAs_fn)
                end


%
%                for k = 1:numel(consess)
%                    str = consess{k}.tcon.name;
%                    montage_fn1 = fullfile(run_dir_stats_new, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-' str '_threshtmap.png']);
%                    montage_fn2 = fullfile(run_dir_stats_new, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-' str '_noFWEthreshtmap.png']);
%                    copyfile(montage_fn1, report_dir);
%                    copyfile(montage_fn2, report_dir)
%                end


%                run_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run '_echo-' echo]);
%                for i = 1:3
%                    fn1 = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', i) '.nii']);
%                    fn2 = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', i) '_binary_clusters.nii']);
%                    fn3 = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', i) '_nary_clusters.nii']);
%                    fnc = fullfile(run_dir_stats, ['con_' sprintf('%04d', i) '.nii']);
%
%                    fnsave1 = strrep(fn1, '.nii', '_MNI152.nii');
%                    fnsave2 = strrep(fn2, '.nii', '_MNI152.nii');
%                    fnsave3 = strrep(fn3, '.nii', '_MNI152.nii');
%                    fnsavec = strrep(fnc, '.nii', '_MNI152.nii');
%
%                    if exist(fn1, 'file')
%                        toTransform_fns = [toTransform_fns, {fn1, fn2, fn3, fnc}];
%                        saveAs_transform_fns = [saveAs_transform_fns, {fnsave1, fnsave2, fnsave3, fnsavec}];
%                    else
%                        continue;
%                    end
%                end

            end
        end
    end

%    transformation_fn = options.indiv_to_mni_fn;
%    template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);
%    fmrwhy_batch_normaliseWrite(toTransform_fns, transformation_fn, template_fn, saveAs_transform_fns)
end
toc;