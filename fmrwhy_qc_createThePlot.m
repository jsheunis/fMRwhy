function fmrwhy_qc_createThePlot(bids_dir, sub, ses, task, run, echo)
%
% This is a Matlab script that uses custom code and spm12 routines to plot
% a version of THE PLOT from Jonathan Power. See:
% - https://www.jonathanpower.net/2017-ni-the-plot.html
% - https://www.ncbi.nlm.nih.gov/pubmed/27510328
% (code has been partly adopted from Power's script)
%
% THE PLOT is a useful visual quality checking tool for fMRI timeseries
% data. It is a 2D plot of voxel intensities over time, with voxels ordered
% in bins (GM, WM, CSF) to indicate some directionality in the data. Motion
% and physiological noise are more easily visually inspected with this
% tool, especially if viewed next to other quality traces, like FD, DVARS
% or physiological recordings.

% 
% UPDATES:
% - 23 July 2019: added voxel ordering options according to RO and GSO
% methods in: https://www.biorxiv.org/content/10.1101/662726v1.full

% Load/create required defaults
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
template_task = 'motor'; % changed for fingertapping experiment. TODO: change back. and update functioning.
template_run = '1';
template_echo = '2';

% Directory and content setup
deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
qc_dir = fullfile(deriv_dir, 'fmrwhy-qc');
sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);
sub_dir_qc = fullfile(qc_dir, ['sub-' sub]);
func_dir_qc = fullfile(sub_dir_qc, 'func');

% Grab anatomical image
anatomical_fn = fullfile(sub_dir_preproc, 'anat', ['sub-' sub '_T1w.nii']);

% Grab template filename
template_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_space-individual_bold.nii']);

% Grab functional timeseries filenames
functional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_bold.nii']);

% Grab basicfunc preproc filenames
motion_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_desc-confounds_motion.tsv']);
afunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-apreproc_bold.nii']);
rfunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rpreproc_bold.nii']);
rafunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_bold.nii']);
sfunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-spreproc_bold.nii']);
srfunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srpreproc_bold.nii']);
srafunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-srapreproc_bold.nii']);
framewise_displacement_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_fd.tsv']);
tissue_regr_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_tissue.tsv']);
physio_regr_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_physio.tsv']);
confounds_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_regressors.tsv']);

% ThePlot settings:
fwhm = 6; % for preprocessing smoothing steps
use_processed = 0;  % use realigned data for the plot? yes = 1
intensity_scale = [-6 6]; % scaling for plot image intensity, see what works
ordering = 1; % ordering of voxels in ThePlotSPM;
% 0 = Random order via standard Matlab indexing (RO)
% 1 = Grey matter signal ordering (GSO) (according to the Pearsons correlation between voxel timeseries and global signal )
% 2 = Cluster-similarity ordering (CO) - NOT IMPLEMENTED YET
% roi = ''; % figure generated for voxels in supplied ROI - NOT IMPLEMENTED YET
if ordering == 0
    order_text = '(RO)';
elseif ordering == 1
    order_text = '(GSO)';
elseif ordering == 2
    order_text = '(CO)';
else
end

% Timeseries to use for ThePlot
functional4D_fn = sfunctional_fn;

% Get image information from
func_spm = spm_vol(functional4D_fn);
tsize = size(func_spm);
Nt = tsize(1);
Ni= func_spm(1).dim(1);
Nj= func_spm(1).dim(2);
Nk= func_spm(1).dim(3);

% Load multiple confound regressors
confounds_struct = tdfread(confounds_fn);
confounds_mat = struct2array(confounds_struct);

% Load mask data
masks = fmrwhy_util_loadMasks(bids_dir, sub);

% Calculate stats data (including BOLD percentage signal change)
stats = fmrwhy_qc_calculateStats(bids_dir, sub, functional4D_fn);

% Load data
F_2D_psc = reshape(stats.data_4D_psc, Ni*Nj*Nk, Nt); %[voxels, time]
I_GM = masks.GM_mask_I;
I_WM = masks.WM_mask_I;
I_CSF = masks.CSF_mask_I;
I_brain = masks.brain_mask_I;

% Voxel ordering
switch(ordering)
    case 0
        GM_img = F_2D_psc(I_GM, :);
        WM_img = F_2D_psc(I_WM, :);
        CSF_img = F_2D_psc(I_CSF, :);
        all_img = [GM_img; WM_img; CSF_img];
    case 1
        GM_signals = F_2D_psc(I_GM, :);
        GM_time_series = mean(GM_signals, 1);
        F_masked_psc = F_2D_psc(I_brain, :);
        Rcorr = corr(F_masked_psc', GM_time_series');
        [R_sorted, I_sorted] = sort(Rcorr,'descend');
        all_img = F_masked_psc(I_sorted, :);
    case 2
        % TODO: implement this correctly
%        % From: https://gist.github.com/benfulcher/f243bdc4ab80a7351f083f35bdd4db63
%        dataMatrix = F_masked_psc';
%        distanceMetric = 'Euclidean';
%        linkageMethod = 'average';
%%         dataMatrixNorm = zscore(dataMatrix); % normalize columns
%        dataMatrixNorm = dataMatrix;
%        R = pdist(dataMatrixNorm,distanceMetric); % Pairwise distances
%        links = linkage(R,linkageMethod); % Do hierarchical linkage
%        f = figure('color','w');
%        set(gcf,'Visible','off'); % suppress figure output
%        [~,~,ord] = dendrogram(links,0);
%        close(f); % close the invisible figure used for the dendrogram
%        dataMatrixClustered = dataMatrixNorm(ord,:); % Reorder rows
%        all_img = dataMatrixClustered';
    otherwise
        GM_img = F_2D_psc(I_GM, :);
        WM_img = F_2D_psc(I_WM, :);
        CSF_img = F_2D_psc(I_CSF, :);
        all_img = [GM_img; WM_img; CSF_img];
end

% Get screen size for plotting
scr_size = get(0,'ScreenSize');
dist = scr_size(4);
if scr_size(3) < dist
    dist = scr_size(3);
end

%xx = [0 0 1 1]; yy1 = [0 30 30 0]; yy2 = [30 60 60 30]; yy3 = [60 100 100 60];
%>> p1 = patch(xx, yy1,'blue','EdgeColor','none')


% Figure
f = figure('units','normalized','outerposition',[0 0 1 1]);
fontsizeL = 15;
fontsizeM = 13;
ax1 = subplot(5,1,[1:4]);
imagesc(ax1, all_img); colormap(gray); caxis(intensity_scale);
title(ax1, ['The Plot SPM ' order_text],'fontsize',fontsizeL)
ylabel(ax1, 'Voxels','fontsize',fontsizeM)
set(ax1,'Xticklabel',[]);
% Add patches if ordered by tissue compartment
if ordering == 0
    hold on;
    gm_pos = numel(I_GM);
    gwm_pos = numel(I_GM) + numel(I_WM);
    gwmcsf_pos = numel(I_GM) + numel(I_WM) + numel(I_CSF);
    x_gm = [0 0 1 1]; y_gm = [0 gm_pos gm_pos 0];
    x_wm = [0 0 1 1]; y_wm = [gm_pos gwm_pos gwm_pos gm_pos];
    x_csf = [0 0 1 1]; y_csf = [gwm_pos gwmcsf_pos gwmcsf_pos gwm_pos];
    p1 = patch(x_gm, y_gm,'green','EdgeColor','none');
    p2 = patch(x_wm, y_wm,'blue','EdgeColor','none');
    p3 = patch(x_csf, y_csf,'red','EdgeColor','none');
    hold off;
end
xlim(ax1,[0 Nt])

ax2 = subplot(5,1,5);
plot(ax2, confounds_struct.framewise_displacement, 'LineWidth', 2); grid;

title(ax2, 'FD','fontsize',fontsizeL)
ylabel(ax2, 'mm','fontsize',fontsizeM)
xlabel(ax2, 'fMRI volumes','fontsize',fontsizeM)

xlim(ax2,[0 Nt])
ylim(ax2,[0 2])
% linkaxes([ax1, ax2], 'x')

theplot_fn = fullfile(func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-RO_grayplot.png']);
if ~exist(theplot_fn, 'file')
    print(f,theplot_fn,'-dpng')
end





