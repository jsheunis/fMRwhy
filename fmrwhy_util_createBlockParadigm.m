function [task_time_course, convolved_ttc] = fmrwhy_util_createBlockParadigm(Nvol, TR, onsets, durations, precision)

% unit = 'scans' / 'seconds'
% onsets = M x 1
% durations = M x 1

%precision = 0.001;
%precision = 1;
task_time_course = zeros(Nvol*TR/precision, 1);

onsets_adj = onsets/precision;
duration_adj = durations/precision;

for i = 1:numel(onsets_adj)
    task_time_course(onsets_adj(i):onsets_adj(i)+duration_adj(i)-1) = 1;
end

% Calculate HRF-convolved task and baseline designs. The HRF (hemodynamic
% response function) tells us that blood (and hence oxygen level) takes
% some seconds to respond to stimulus-driven neuronal activity. Thus when
% we e.g. see a stimulus, the oxygen level change due to the neuronal
% activity will only be measurable with fMRI a few seconds later. So we do
% convolution of the task design to have a more accurate expectatioon for
% the fMRI signal change.
hrf = spm_hrf(TR);
convolved_ttc = spm_Volterra(struct('u', task_time_course, 'name', {{'task'}}), hrf);
convolved_ttc = convolved_ttc(1:Nvol);


% Condition-encoded task-baseline vector
%vectEncCond = task_design + 1; % 1 = index in baseline block; 2 = index in task block
%% Get indices for task and baseline blocks
%task_blocks = bwlabel(task_design);
%baseline_blocks = bwlabel(~task_design);
%baseline_design = (~task_design);
%I_baseline = [];
%for n = 1:length(baseline_onsets)
%    I_block = find(baseline_blocks == n);
%    if n > 1
%        I_block = I_block((1+TR_skip):end);
%    end
%    I_baseline = [I_baseline I_block];
%end