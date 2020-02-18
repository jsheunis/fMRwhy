function PSC = getPercentageSignalChange(image, N_task_blocks, N_rest_blocks, block_duration, task_onsets, Nt)
% PSC = getPercentageSignalChange(6, 7, 16, [17;49;81;113;145;177], 208)


% ...
% this function needs to be much better
% ...


% gets PSC per volume (only using task volumes) for 4D image. PSC from mean
% of resting block preceding task

% N_task_blocks = 6;
% N_rest_blocks = 7;
% block_duration = 16; % number of volumes per task/rest block
% task_onsets = [17;49;81;113;145;177]; % starting volumes for task blocks

task_design = zeros(1,Nt);

for i = 1:numel(task_onsets)
    for j = 1:block_duration
        task_design(1,task_onsets(i)+j-1) = 1;
    end
end

rest_design = ~task_design;
rest_blocks = bwlabel(rest_design);
task_blocks = bwlabel(task_design);
rest_mean = cell(N_rest_blocks,2);

% calculate rest means
for i = 1:N_rest_blocks
    
    rest_mean{i,1} = mean(image(:,find(rest_blocks==i)),2);
    if i ~= N_rest_blocks
        rest_mean{i,2} = repmat(rest_mean{i,1}, 1, 2*block_duration);
    else
        rest_mean{i,2} = repmat(rest_mean{i,1}, 1, block_duration);
    end
    
    if i==1
        mean_blocks = rest_mean{i,2};
    else
        mean_blocks = cat(2, mean_blocks, rest_mean{i,2});
    end
end
PSC = struct;
PSC.PSC_blocks = 100*(image./mean_blocks) - 100;
PSC.PSC_blocks(isnan(PSC.PSC_blocks))=0;

PSC.mean_blocks = mean_blocks;
PSC.rest_mean = rest_mean;
PSC.rest_blocks = rest_blocks;
PSC.task_blocks = task_blocks;
PSC.task_design = task_design;



