% Generic run-through

data_dir = '/Users/jheunis/Documents/MATLAB/ME/Volunteer_data/V10all';
data_dir_new = '/Users/jheunis/Documents/MATLAB/ME/Volunteer_data/V10all';
subj_list = [1 2];
subj_runs = [7 8];
info_mat = zeros(75, 13);
% subj run pipeline T-threshold Nclusters Nvoxels Tmax Tmax-voxel-location*3
% Tmax-mm-location*3
j = 0;

for subj = 1:numel(subj_list)
    for run = 1:subj_runs(subj)
        subj_dir = ['subj_' num2str(subj_list(subj)) '_' num2str(run)];
        disp(subj_dir);
        N_pipelines = 5;
        pipeline_dir = cell(N_pipelines,1);
        pipeline_dir{1} = [data_dir filesep subj_dir filesep 'Task' filesep 'Pipeline 1 - Kundu'];
        pipeline_dir{2} = [data_dir filesep subj_dir filesep 'Task' filesep 'Pipeline 2 - TE2'];
        pipeline_dir{3} = [data_dir filesep subj_dir filesep 'Task' filesep 'Pipeline 3 - Combined_TE'];
        pipeline_dir{4} = [data_dir filesep subj_dir filesep 'Task' filesep 'Pipeline 4 - Flattened_S0'];
        pipeline_dir{5} = [data_dir filesep subj_dir filesep 'Task' filesep 'Pipeline 5 - T2star'];
        mni_stats_dir = cell(N_pipelines,1);
        main_string_size = numel([data_dir_new filesep subj_dir filesep 'Task' filesep 'Pipeline X - ']);
        
        for i = 1:N_pipelines
            j = j+1;
            mni_stats_dir{i} = [pipeline_dir{i} filesep 'MNIstats'];
            load([mni_stats_dir{i} filesep 'xSPM.mat'], 'xSPM')
            pause(1);
            height_threshold = xSPM.u;
            [Tmax, index] = max(xSPM.Z);
            [N, Nvoxels] = size(xSPM.Z);
            Vx= xSPM.XYZ(1,index);
            Vy = xSPM.XYZ(2,index);
            Vz = xSPM.XYZ(3,index);
            Vxmm= xSPM.XYZmm(1,index);
            Vymm = xSPM.XYZmm(2,index);
            Vzmm = xSPM.XYZmm(3,index);
            tmap = spm_read_vols(spm_vol([mni_stats_dir{i} filesep 'spmT_0001.nii']));
            tmap_thresholded = tmap;
            tmap_thresholded(tmap_thresholded<height_threshold) = 0;
            tmap_thresholded(tmap_thresholded>0) = 1;
            [L,Nclusters] = spm_bwlabel(tmap_thresholded,18);
            info_mat(j, :) = [subj run i height_threshold Nclusters Nvoxels Tmax Vx Vy Vz Vxmm Vymm Vzmm];
        end
    end
end