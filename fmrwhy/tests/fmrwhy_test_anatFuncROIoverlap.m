bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';
sub = '001';
ses = '';
task = 'motor';
run = '2';
echo = '2';

options = struct;

options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% Update workflow params with subject functional derivative filenames
options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, options.template_echo, options);

roi_fn = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-rleftMotor_roi.nii']);
% anat_roi_fn{1} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-r' options.roi.(task).desc{1} '_roi.nii']);
fn1 = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-motor_run-1_echo-2_noFWEp001e20/spmT_0001_binary_clusters.nii';

output = fmrwhy_util_getBinaryOverlap({roi_fn, fn1});

%
% run_dir_stats = fullfile(options.sub_dir_stats, ['task-' task '_run-' run]);
% consess = options.firstlevel.(task).(['run' run]).contrast_params.consess;
% background_fn = options.rcoregest_anatomical_fn;
% [p, frm, rg, dim] = fmrwhy_util_readOrientNifti(background_fn);
% background_img = p.nii.img;
% anat_roi_fn = {};
% anat_roi_fn{1} = fullfile(options.anat_dir_preproc, ['sub-' sub '_space-individual_desc-r' options.roi.(task).desc{1} '_roi.nii']);
%
%% For each anatomical ROI:
% for j = 1:numel(anat_roi_fn)
%    % Grab binary image of ROI, not yet for plotting
%    nii = nii_tool('load', anat_roi_fn{j});
%    anat_roi_img = fmrwhy_util_createBinaryImg(nii.img, 0.1);
%    % For each task contrast, i.e. for each tmap, grab binary tmap clusters, find overlap with anatomical ROI, save
%    for k = 1:numel(consess)
%%        tmap_fn = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', k) '.nii']);
%        tmap_clusters_fn = fullfile(run_dir_stats, ['spmT_' sprintf('%04d', k) '_binary_clusters.nii']);
%%        [ptmap, ~, ~, ~] = fmrwhy_util_readOrientNifti(tmap_fn);
%        nii = nii_tool('load', tmap_clusters_fn);
%        func_roi_img = fmrwhy_util_createBinaryImg(nii.img, 0.1);
%
%        overlap_img = anat_roi_img & func_roi_img;
%        size_anat_roi = sum(anat_roi_img(:))
%        size_func_roi = sum(func_roi_img(:))
%        size_overlap = sum(overlap_img(:))
%        dice_coeff = 2*size_overlap/(size_anat_roi + size_func_roi)
%
%        no_scaling = 1;
%        str = consess{k}.tcon.name;
%        saveAs_overlap_fn = fullfile(run_dir_stats, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-' str 'Overlaps' options.roi.(task).desc{j} '_roi.nii']);
%        fmrwhy_util_saveNifti(saveAs_overlap_fn, double(overlap_img), tmap_clusters_fn, no_scaling);
%        [poverlap, ~, ~, ~] = fmrwhy_util_readOrientNifti(saveAs_overlap_fn);
%        plot_overlap_img = double(poverlap.nii.img);
%        saveAs_fn = strrep(saveAs_overlap_fn, '.nii', '.png');
%        roi_rgbcolors = [255, 0, 0];
%
%%        [ptmapc, ~, ~, ~] = fmrwhy_util_readOrientNifti(tmap_clusters_fn);
%%        func_roi_img = double(ptmapc.nii.img);
%%
%%        [panat, ~, ~, ~] = fmrwhy_util_readOrientNifti(anat_roi_fn{j});
%%        anat_roi_img = fmrwhy_util_createBinaryImg(panat.nii.img, 0.1);
%
%%        funcmontage = fmrwhy_util_createStatsOverlayMontage(background_img, [], func_roi_img, 9, 1, '', 'gray', 'on', 'max', [], [], roi_rgbcolors, false, [])
%%        anatmontage = fmrwhy_util_createStatsOverlayMontage(background_img, [], anat_roi_img, 9, 1, '', 'gray', 'on', 'max', [], [], roi_rgbcolors, false, [])
%%        overlapss = fmrwhy_util_createStatsOverlayMontage(background_img, [], double(overlap_img), 9, 1, '', 'gray', 'on', 'max', [], [], roi_rgbcolors, false, [])
%        overlapmontage = fmrwhy_util_createStatsOverlayMontage(background_img, [], plot_overlap_img, 9, 1, '', 'gray', 'off', 'max', [], [], roi_rgbcolors, false, saveAs_fn);
%
%    end
% end
