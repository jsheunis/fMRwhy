function output = fmrwhy_qc_createMaskMontages(bids_dir, sub, savefigs)
% Function to create montages of tissue mask contours overlaid on
% assumes standard fmrwhy-preproc directory structure
% assumes fmrwhy_preproc_structFunc.m has been run successfully

% Template values
template_task = 'rest'; % TODO: updated for fingertapping experiment, change!!!
template_run = 1;

% BIDS structure values
BIDS = spm_BIDS(bids_dir);
subjects = spm_BIDS(BIDS,'subjects');
sessions = spm_BIDS(BIDS,'sessions');
runs = spm_BIDS(BIDS,'runs');
tasks = spm_BIDS(BIDS,'tasks');
types = spm_BIDS(BIDS,'types');
modalities = spm_BIDS(BIDS,'modalities');

% Directory and content setup
deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
qc_dir = fullfile(deriv_dir, 'fmrwhy-qc');
sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);
sub_dir_qc = fullfile(qc_dir, ['sub-' sub]);
anat_dir_qc = fullfile(sub_dir_qc, 'anat');
if ~exist(anat_dir_qc, 'dir')
    mkdir(anat_dir_qc)
end

% Get functional template
template_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' num2str(template_run) '_space-individual_bold.nii']);
template_img = spm_read_vols(spm_vol(template_fn));

% Get anatomical masks in individual functional space
masks = fmrwhy_util_loadMasks(bids_dir, sub);
mask_images = {masks.GM_mask_3D, masks.WM_mask_3D, masks.CSF_mask_3D, masks.brain_mask_3D};
mask_names = {'Grey matter', 'White matter', 'Cerebrospinal fluid', 'Whole brain'};
mask_montage_fns = {'_GM_mask_montage', '_WM_mask_montage', '_CSF_mask_montage', '_brain_mask_montage'};
for i = 1:numel(mask_montage_fns)
    mask_montage_fns{i} = fullfile(anat_dir_qc, ['sub-' sub mask_montage_fns{i} '.png']);
end

% Structure to save output
output = struct;

% Create background montage
montage_EPI = fmrwhy_util_createMontage(template_img, 9, 1, 'Template volume', 'gray', 'off', 'max');

% Get screen size for plotting
scr_size = get(0,'ScreenSize');
dist = scr_size(4);
if scr_size(3) < dist
    dist = scr_size(3);
end

% Create figures with background montage and overlaid masks
% TODO: create a sub function for this, i.e. overlaying montages
for i = 1:numel(mask_images)
    %if ~exist(mask_montage_fns{i}, 'file')
        %f(i) = figure('units','pixels','outerposition',[0 0 dist dist]);
        f(i) = figure('units','normalized','outerposition',[0 0 1 1]);
        im1 = imagesc(montage_EPI.whole_img);
        colormap('gray');
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset; 
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
        hold(ax, 'on')
        [Nimx, Nimy] = size(montage_EPI.whole_img);
        oo = ones(Nimx, Nimy);
        zz = zeros(Nimx, Nimy);
        red = cat(3, oo, zz, zz);
        montage_mask = fmrwhy_util_createMontage(mask_images{i}, 9, 1, mask_names{i}, 'gray', 'off', 'max');
        bound_whole_bin = bwboundaries(montage_mask.whole_img);
        Nblobs_bin = numel(bound_whole_bin);
        for b = 1:Nblobs_bin
            p = plot(ax, bound_whole_bin{b,1}(:,2), bound_whole_bin{b,1}(:,1), 'r', 'LineWidth', 1);
        end
        imC = imagesc(ax, red);
        set(imC, 'AlphaData', 0.2*montage_mask.whole_img);
        hold(ax, 'off');
        set(ax,'xtick',[])
        set(ax,'xticklabel',[])
        set(ax,'ytick',[])
        set(ax,'yticklabel',[])
        set(ax,'ztick',[])
        set(ax,'zticklabel',[])
        print(f(i),mask_montage_fns{i},'-dpng', '-r0')
    %end
end




% TODO: incorporate a better way to show/hide figures and to save images or not, using function arguments

