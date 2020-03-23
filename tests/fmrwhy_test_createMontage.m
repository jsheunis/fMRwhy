

template_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/func/sub-001_task-rest_run-1_space-individual_bold.nii';
anat_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_T1w.nii';
overlay_fn = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-preproc/sub-001/anat/sub-001_space-individual_desc-rrightMotor_roi.nii';


[p1, frm1, rg1, dim1] = fmrwhy_util_readNifti(template_fn);
[p2, frm2, rg2, dim2] = fmrwhy_util_readNifti(anat_fn);
[p3, frm3, rg3, dim3] = fmrwhy_util_readNifti(overlay_fn);

overlay_img = fmrwhy_util_createBinaryImg(p3.nii.img, 0); % left_motor_roi_img

rotate = 1;
montage1 = fmrwhy_util_createMontage(p1.nii.img, 9, rotate, 'test', 'gray', 'off', 'max');
montage2 = fmrwhy_util_createMontage(p2.nii.img(:,:,181:214), 9, rotate, 'test', 'gray', 'off', 'max');


% Create figures with background montage and overlaid masks
%f(i) = figure('units','pixels','outerposition',[0 0 dist dist]);
f = figure('units','normalized','outerposition',[0 0 1 1]);
im1 = imagesc(montage1.whole_img);
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
[Nimx, Nimy] = size(montage1.whole_img);
oo = ones(Nimx, Nimy);
zz = zeros(Nimx, Nimy);
red = cat(3, oo, zz, zz);
montage_overlay = fmrwhy_util_createMontage(overlay_img, 9, rotate, 'Overlay', 'gray', 'off', 'max');
bound_whole_bin = bwboundaries(montage_overlay.whole_img);
Nblobs_bin = numel(bound_whole_bin);
for b = 1:Nblobs_bin
p = plot(ax, bound_whole_bin{b,1}(:,2), bound_whole_bin{b,1}(:,1), 'r', 'LineWidth', 1);
end
imC = imagesc(ax, red);
set(imC, 'AlphaData', 0.2*montage_overlay.whole_img);
hold(ax, 'off');
%set(ax,'xtick',[])
%set(ax,'xticklabel',[])
%set(ax,'ytick',[])
%set(ax,'yticklabel',[])
%set(ax,'ztick',[])
%set(ax,'zticklabel',[])


f = figure('units','normalized','outerposition',[0 0 1 1]);
im2 = imagesc(montage2.whole_img);
colormap('gray');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

