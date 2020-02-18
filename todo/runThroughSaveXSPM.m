% Generic run-through

data_dir = '/Users/jheunis/Documents/MATLAB/ME/Volunteer_data/V10all';
subj_list = [1 2];
subj_runs = [7 8];

for subj = 2:2%numel(subj_list)
    for run = 5:5%subj_runs(subj)
        subj_dir = ['subj_' num2str(subj_list(subj)) '_' num2str(run)];
        
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
            if i == 1
                png_fn{i} = [data_dir filesep subj_dir '_TE2_orig.png'];
            else
                png_fn{i} = [data_dir filesep subj_dir '_' pipeline_dir{i}((main_string_size+1):end) '.png'];
            end
            
            mni_stats_dir{i} = [pipeline_dir{i} filesep 'MNIstats'];
            viewSpmStatsResults(jobs_dir, mni_stats_dir{i})
            pause(5);
            save([mni_stats_dir{i} filesep 'xSPM.mat'], 'xSPM')
            pause(1);
        end
    end
end