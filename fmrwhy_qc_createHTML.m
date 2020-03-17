function output = fmrwhy_qc_createHTML(bids_dir, sub, ses, task, run)

% Load/create required defaults
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
template_task = 'motor'; % changed for fingertapping experiment. TODO: change back. and update functioning.
template_run = '1';
template_echo = '2';
stre_txt = ['sub-' sub '_task-' task '_run-' run '_echo-' echo];
str_txt = ['sub-' sub '_task-' task '_run-' run];

% Directory and content setup
deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
qc_dir = fullfile(deriv_dir, 'fmrwhy-qc');
sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);
sub_dir_qc = fullfile(qc_dir, ['sub-' sub]);
func_dir_qc = fullfile(sub_dir_qc, 'func');
anat_dir_qc = fullfile(sub_dir_qc, 'anat');

% Create HTML log file
dt = datetime('now');
[Y,MO,D,H,MI,S] = datevec(dt);
dt_str = [num2str(Y) num2str(MO) num2str(D) num2str(H) num2str(MI) num2str(round(S))];
t = datestr(dt);

log_name = [subject '_' dt_str '.html'];

fid = fopen(log_name,'a');
fprintf(fid, '<H2>Log</H2>');
fprintf(fid, ['\n<BR>Subject:  ' subject]);

fprintf(fid, ['\n<BR>Date/time:  ' t]);

fprintf(fid, '<H2>Imaging info</H2>');
fprintf(fid, ['\nVolumes:  ' num2str(Nt)]);
fprintf(fid, ['\n<BR>Voxels (x,y,z):  ' num2str(Ni) ', ' num2str(Nj) ', ' num2str(Nk)]);

fprintf(fid, '<H2>Timeseries summary</H2>');
fprintf(fid, '\n<TABLE><TR><TD><img src="timeseries_summary.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );

fprintf(fid, '<H2>QC Metrics</H2>');
fprintf(fid, ['\nFD threshold (mm):  ' num2str(FD_threshold)]);
fprintf(fid, ['\n<BR>FD outliers:  ' num2str(numel(FD_measures.FD_outliers_ind))]);
fprintf(fid, ['\n<BR>Total FD:  ' num2str(FD_measures.FD_sum)]);
fprintf(fid, ['\n<BR>Mean FD:  ' num2str(FD_measures.FD_mean)]);
fprintf(fid, ['\n<BR>Mean Zscore:  ' num2str(Zstat_mean)]);
fprintf(fid, ['\n<BR>GCOR:  ' num2str(GCOR)]);
fprintf(fid, ['\n<BR>tSNR (brain):  ' num2str(tSNR_brain)]);
fprintf(fid, ['\n<BR>tSNR (GM):  ' num2str(tSNR_GM)]);
fprintf(fid, ['\n<BR>tSNR (WM):  ' num2str(tSNR_WM)]);
fprintf(fid, ['\n<BR>tSNR (CSF):  ' num2str(tSNR_CSF)]);

fprintf(fid, '<H2>Anatomical QC</H2>');
fprintf(fid, '\n<TABLE><TR><TD><img src="mean_epi.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );
fprintf(fid, '\n<TABLE><TR><TD><img src="stddev_epi.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );
fprintf(fid, '\n<TABLE><TR><TD><img src="tsnr_whole.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );
fprintf(fid, '\n<TABLE><TR><TD><img src="tsnr_brain.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );
fprintf(fid, '\n<TABLE><TR><TD><img src="mask_contour.png" alt="no picture" width=700 height=700></TD></TR></TABLE>' );
fclose(fid);


























