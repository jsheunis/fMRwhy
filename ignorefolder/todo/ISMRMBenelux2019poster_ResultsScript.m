% SCRIPT: ISMRMBenelux2019poster_ResultsScript
%--------------------------------------------------------------------------
% Copyright (C) Neu3CA Research Group, Eindhoven University of Technology
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,USA.
%
% Author: Stephan Heunis, <j.s.heunis@tue.nl>, 2019

%--------------------------------------------------------------------------
% DEFINITION
%--------------------------------------------------------------------------

% The "ISMRMBenelux2019poster_ResultsScript.m" is called by the script:
% "ISMRMBenelux2019poster_ProcessingScript.m". It uses all data generated
% by the previous script in order to generate the results shown in the
% poster presented at the 11th Annual meeting of the ISMRM Benelux chapter.
% See: https://doi.org/10.5281/zenodo.2553256

% The script requires the following setup/installs:
%   -   Successful execution of ISMRMBenelux2019poster_ProcessingScript.m,
%       with all variables still active in the Matlab workspace.
%   -   Matlab 2014a or later
%   -   SPM12: http://www.fil.ion.ucl.ac.uk/spm/
%   -   Raincloudplots Matlab toolbox:
%       https://github.com/RainCloudPlots/RainCloudPlots (with its own dependencies)

%% Initialize variables

MNI_size = [79 95 79];
Nvox_mni = MNI_size(1)*MNI_size(2)*MNI_size(3);
results = cell(7,3);
for i = 1:7
    results{i,1} = zeros(Nvox_mni, Nsub); % whole image
    results{i,2} = zeros(Nvox_mni, Nsub); % whole brain
    results{i,3} = zeros(Nvox_mni, Nsub); % grey matter
end
subj_summary = zeros(Nsub, 14);
tsnr_threshold = 250;
perc_threshold = 250;
all_subs_brain_mask = false(Nvox_mni, Nsub);
all_subs_GM_mask = false(Nvox_mni, Nsub);

%% For each subject:
%   load processed data
%   get masks in MNI space
%   calculate tSNR percentage difference maps
%   calculate subject summary values

for subj = 1:Nsub
    
    % Initialize subject file names
    subj_dir = ['sub-' sprintf('%02d', subj)];
    disp(subj_dir);
    anat_dir = [data_dir filesep subj_dir filesep 'anat'];
    func_dir = [data_dir filesep subj_dir filesep 'func'];
    s_fn = [anat_dir filesep subj_dir '_T1W.nii'];
    cd([data_dir filesep subj_dir]);
    
    % LOAD PREPROCESSED data
    disp('Load preprocessed data')
    % Preprocess structural and f0 images
    [d, f, e] = fileparts(s_fn);
    preproc_data = struct;
    preproc_data.forward_transformation = [d filesep 'y_' f e];
    preproc_data.inverse_transformation = [d filesep 'iy_' f e];
    preproc_data.gm_fn = [d filesep 'c1' f e];
    preproc_data.wm_fn = [d filesep 'c2' f e];
    preproc_data.csf_fn = [d filesep 'c3' f e];
    preproc_data.bone_fn = [d filesep 'c4' f e];
    preproc_data.soft_fn = [d filesep 'c5' f e];
    preproc_data.air_fn = [d filesep 'c6' f e];
    preproc_data.rstructural_fn = [d filesep 'r' f e];
    preproc_data.rgm_fn = [d filesep 'rc1' f e];
    preproc_data.rwm_fn = [d filesep 'rc2' f e];
    preproc_data.rcsf_fn = [d filesep 'rc3' f e];
    preproc_data.rbone_fn = [d filesep 'rc4' f e];
    preproc_data.rsoft_fn = [d filesep 'rc5' f e];
    preproc_data.rair_fn = [d filesep 'rc6' f e];
    preproc_data.wrstructural_fn = [d filesep 'wr' f e];
    preproc_data.wrgm_fn = [d filesep 'wrc1' f e];
    preproc_data.wrwm_fn = [d filesep 'wrc2' f e];
    preproc_data.wrcsf_fn = [d filesep 'wrc3' f e];
    preproc_data.wrbone_fn = [d filesep 'wrc4' f e];
    preproc_data.wrsoft_fn = [d filesep 'wrc5' f e];
    preproc_data.wrair_fn = [d filesep 'wrc6' f e];
    
    % Construct MNI space GM, WM and CSF masks (separately and together)
    [GM_img_bin, WM_img_bin, CSF_img_bin] = createBinarySegments(preproc_data.wrgm_fn, preproc_data.wrwm_fn, preproc_data.wrcsf_fn, 0.5);
    I_GM = find(GM_img_bin);
    I_WM = find(WM_img_bin);
    I_CSF = find(CSF_img_bin);
    I_total_test = numel(I_GM)+numel(I_WM)+numel(I_CSF);
    mask_reshaped = GM_img_bin | WM_img_bin | CSF_img_bin;
    I_mask = find(mask_reshaped);
    Nmaskvox = numel(I_mask);
    Nvox = numel(GM_img_bin);
    [Ni, Nj, Nk] = size(GM_img_bin);
    
    all_subs_brain_mask(:, subj) = mask_reshaped(:);
    all_subs_GM_mask(:, subj) = GM_img_bin(:);
    
    % LOAD TSNR MAPS
    disp('Load tSNR maps')
    % 1 - 2nd echo realigned
    % 2 - combined echoes using precalculated T2star for weighting
    % 3 - combined echoes using precalculated tSNR for weighting
    % 4 - combined echoes using real-time T2star for weighting
    tSNR = cell(4,1);
    tSNR_img = cell(4,1);
    for i = 1:numel(tSNR)
        tsnr_fn = [data_dir filesep subj_dir filesep 'wtSNR_ts' num2str(i) '_sub-' sprintf('%02d', subj) '.nii'];
        tSNR_img{i} = spm_read_vols(spm_vol(tsnr_fn));
        tSNR{i} = reshape(tSNR_img{i}, Ni*Nj*Nk, 1);
    end

    % Calculate percentage difference maps for j=1:3 (whole image, whole brain mask, grey matter mask)
    disp('Calculate percentage difference maps')
    tSNR_diffs =  cell(3,1);
    tSNR_diffs_img =  cell(3,1);
    tSNR_perc =  cell(3,1);
    tSNR_perc_img =  cell(3,1);
    mask = 0;
    j = 1;
    for i = 1:3
        tSNR_diffs{i} = zeros(Nvox,1);
        tSNR_diffs{i} = tSNR{i+1} - tSNR{1};
        tSNR_perc{i} = zeros(Nvox,1);
        tSNR_perc{i} = tSNR_diffs{i}./tSNR{1}*100;
        tSNR_diffs_img{i} = reshape(tSNR_diffs{i}, Ni, Nj, Nk);
        tSNR_perc_img{i} = reshape(tSNR_perc{i}, Ni, Nj, Nk);
    end
    results{1,j}(:,subj) =  tSNR{1};
    results{2,j}(:,subj) =  tSNR{2};
    results{3,j}(:,subj) =  tSNR{3};
    results{4,j}(:,subj) =  tSNR{4};
    results{5,j}(:,subj) =  tSNR_perc{1};
    results{6,j}(:,subj) =  tSNR_perc{2};
    results{7,j}(:,subj) =  tSNR_perc{3};
    
    mask = I_mask;
    j = 2;
    for i = 1:3
        tSNR_diffs{i} = zeros(Nvox,1);
        tSNR_diffs{i}(mask) = tSNR{i+1}(mask) - tSNR{1}(mask);
        tSNR_perc{i} = zeros(Nvox,1);
        tSNR_perc{i}(mask) = tSNR_diffs{i}(mask)./tSNR{1}(mask)*100;
        tSNR_diffs_img{i} = reshape(tSNR_diffs{i}, Ni, Nj, Nk);
        tSNR_perc_img{i} = reshape(tSNR_perc{i}, Ni, Nj, Nk);
    end
    results{1,j}(:,subj) =  tSNR{1};
    results{2,j}(:,subj) =  tSNR{2};
    results{3,j}(:,subj) =  tSNR{3};
    results{4,j}(:,subj) =  tSNR{4};
    results{5,j}(:,subj) =  tSNR_perc{1};
    results{6,j}(:,subj) =  tSNR_perc{2};
    results{7,j}(:,subj) =  tSNR_perc{3};
    
    mask = I_GM;
    j = 3;
    for i = 1:3
        tSNR_diffs{i} = zeros(Nvox,1);
        tSNR_diffs{i}(mask) = tSNR{i+1}(mask) - tSNR{1}(mask);
        tSNR_perc{i} = zeros(Nvox,1);
        tSNR_perc{i}(mask) = tSNR_diffs{i}(mask)./tSNR{1}(mask)*100;
        tSNR_diffs_img{i} = reshape(tSNR_diffs{i}, Ni, Nj, Nk);
        tSNR_perc_img{i} = reshape(tSNR_perc{i}, Ni, Nj, Nk);
    end
    results{1,j}(:,subj) =  tSNR{1};
    results{2,j}(:,subj) =  tSNR{2};
    results{3,j}(:,subj) =  tSNR{3};
    results{4,j}(:,subj) =  tSNR{4};
    results{5,j}(:,subj) =  tSNR_perc{1};
    results{6,j}(:,subj) =  tSNR_perc{2};
    results{7,j}(:,subj) =  tSNR_perc{3};
    
    mask = I_mask;
    j = 2;
    for i = 1:7
        new_mask = find(isfinite(results{i,j}(mask,subj)));
        subj_summary(subj, i) = mean(results{i,j}(new_mask,subj));
    end
    
    mask = I_GM;
    j = 3;
    for i = 8:14
        new_mask = find(isfinite(results{(i-7),j}(mask,subj)));
        subj_summary(subj, i) = mean(results{(i-7),j}(new_mask,subj));
    end
end

%% Now that all tSNR maps and tSNR percentage difference maps are available for all subjects in MNI space:
% calculate averages of the gorup:

results_means = cell(7,3);
thresh = 10000;
for i = 1:7
   results_means{i,1} = thresholdVector(mean(results{i,1}, 2), 0, thresh, 0);
   results_means{i,2} = thresholdVector(mean(results{i,2}, 2), 0, thresh, 0);
   results_means{i,3} = thresholdVector(mean(results{i,3}, 2), 0, thresh, 0);
end

all_GM = find(all(all_subs_GM_mask, 2));
all_brain = find(all(all_subs_brain_mask, 2));

%% Display average tSNR maps and percentage difference maps (in montage format)
tsnr_colormap = 'hot';
perc_colormap = 'parula';
mask_type = 1;
j = 2;

results_means_ave = cell(7,1);

for i = 1:7
    results_means_ave{i} = zeros(592895, 1);
    results_means_ave{i}(all_brain) = results_means{i,j}(all_brain);
end

% For Figure 3 A-D in poster
montage1 = createMontage(reshape(results_means_ave{1}, Ni, Nj, Nk), 8, 1, 'Average tSNR - Echo 2', tsnr_colormap, 0, [0 tsnr_threshold]);
montage2 = createMontage(reshape(results_means_ave{2}, Ni, Nj, Nk), 8, 1, 'Average tSNR - Combined preT2star', tsnr_colormap, 0, [0 tsnr_threshold]);
montage3 = createMontage(reshape(results_means_ave{3}, Ni, Nj, Nk), 8, 1, 'Average tSNR - Combined preTSNR', tsnr_colormap, 0, [0 tsnr_threshold]);
montage4 = createMontage(reshape(results_means_ave{4}, Ni, Nj, Nk), 8, 1, 'Average tSNR - Combined rtT2star', tsnr_colormap, 0, [0 tsnr_threshold]);

% For Figure 4 A-C in poster
montage5 = createMontage(reshape(results_means_ave{5}, Ni, Nj, Nk), 8, 1, 'Average percentage increase - Combined preT2star vs Echo 2', perc_colormap, 0, [0 perc_threshold]);
montage6 = createMontage(reshape(results_means_ave{6}, Ni, Nj, Nk), 8, 1, 'Average percentage increase - Combined preTSNR vs Echo 2', perc_colormap, 0, [0 perc_threshold]);
montage7 = createMontage(reshape(results_means_ave{7}, Ni, Nj, Nk), 8, 1, 'Average percentage increase - Combined rtT2star vs Echo 2', perc_colormap, 0, [0 perc_threshold]);


%% Display raincloudplots - two figures in poster
[cb] = cbrewer('qual','Set3',12,'pchip');
cl(1,:) = cb(4,:);
cl(2,:) = cb(1,:);
fig_position = [100 65 1089 640]; % coordinates for figures

% Figure 5 A and B in poster
f7 = figure('Position', fig_position);
subplot(1,2,1)
set(gca,'XLim', [0 220], 'YLim', [-.03 .04]);
h1= raincloud_plot('X', results_means{1,2}(all_brain) , 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
h2= raincloud_plot('X', results_means{2,2}(all_brain), 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
h3= raincloud_plot('X', results_means{3,2}(all_brain), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
h4= raincloud_plot('X', results_means{4,2}(all_brain), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0);
legend([h1{1} h2{1} h3{1} h4{1}], {'Echo 2', 'pre-T2{\ast}', 'pre-tSNR', 'rt-T2{\ast}'})
title('tSNR distributions within whole brain mask')
set(gca,'XLim', [0 220], 'YLim', [-.03 .04]);
box off
subplot(1,2,2)
set(gca,'XLim', [0 220], 'YLim', [-.03 .04]);
h11= raincloud_plot('X', results_means{1,3}(all_GM) , 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
h22= raincloud_plot('X', results_means{2,3}(all_GM), 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
h33= raincloud_plot('X', results_means{3,3}(all_GM), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
h44= raincloud_plot('X', results_means{4,3}(all_GM), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0);
legend([h11{1} h22{1} h33{1} h44{1}], {'Echo 2', 'pre-T2{\ast}', 'pre-tSNR', 'rt-T2{\ast}'})
title('tSNR distributions within grey matter mask')
set(gca,'XLim', [0 220], 'YLim', [-.03 .04]);
box off

% Figure 6 A and B in poster
f9 = figure('Position', fig_position);
subplot(1,2,1)
h5= raincloud_plot('X', results_means{5,2}(all_brain), 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
h6= raincloud_plot('X', results_means{6,2}(all_brain), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
h7= raincloud_plot('X', results_means{7,2}(all_brain), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
legend([h5{1} h6{1} h7{1}], {'pre-T2{\ast}', 'pre-tSNR', 'rt-T2{\ast}'})
title('Percentage tSNR increase distribution within whole brain mask')
% set(gca,'XLim', [0 240]);
set(gca,'XLim', [0 300], 'YLim', [-.02 .03]);
box off
subplot(1,2,2)
h55 = raincloud_plot('X', results_means{5,3}(all_GM), 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
h66 = raincloud_plot('X', results_means{6,3}(all_GM), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
h77 = raincloud_plot('X', results_means{7,3}(all_GM), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
legend([h55{1} h66{1} h77{1}], {'pre-T2{\ast}', 'pre-tSNR', 'rt-T2{\ast}'})
title('Percentage tSNR increase distribution within grey matter mask')
set(gca,'XLim', [0 300], 'YLim', [-.02 .03]);
box off


%% Extra figure (not in poster)
f10 = figure('Position', fig_position);
h8= raincloud_plot('Y', subj_summary(:,8), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);
h9= raincloud_plot('Y', subj_summary(:,9), 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
h10= raincloud_plot('Y', subj_summary(:,10), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
h11= raincloud_plot('Y', subj_summary(:,11), 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0);
legend([h8{1} h9{1} h10{1} h11{1}], {'Echo 2', 'pre-T2{\ast}', 'pre-tSNR', 'rt-T2{\ast}'})
title('Distributions of average tSNR within whole brain mask')
set(h8{2},'MarkerEdgeColor', 'k', 'SizeData', 25) % handles 1-6 are the cloud area, scatterpoints, and boxplot elements respectively
set(h9{2},'MarkerEdgeColor', 'k', 'SizeData', 25)
set(h10{2},'MarkerEdgeColor', 'k', 'SizeData', 25)
set(h11{2},'MarkerEdgeColor', 'k', 'SizeData', 25)
% set(h2{2}, 'MarkerEdgeColor', 'red') % 
% set(gca,'XLim', [0 240]);
set(gca,'XLim', [0 40], 'YLim', [-.15 .15]);
view([90 -90])
box off
