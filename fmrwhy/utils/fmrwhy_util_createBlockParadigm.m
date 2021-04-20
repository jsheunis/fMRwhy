function [task_time_course, convolved_ttc] = fmrwhy_util_createBlockParadigm(Nvol, TR, onsets, durations, precision, timing_units)

    % unit = 'scans' / 'secs'
    % onsets = M x 1
    % durations = M x 1

    % precision = 0.001;
    % precision = 1;

    % The HRF (hemodynamic response function) tells us that blood (and hence oxygen level)
    % takes some seconds to respond to stimulus-driven neuronal activity. Thus when
    % we e.g. see a stimulus, the oxygen level change due to the neuronal
    % activity will only be measurable with fMRI a few seconds later. So we do
    % convolution of the task design to have a more accurate expectation for
    % the fMRI signal change.

    if nargin < 6
        timing_units = 'scans';
    end

    if strcmp(timing_units, 'scans')
        task_time_course = zeros(Nvol * TR / precision, 1);
        onsets_adj = onsets / precision;
        duration_adj = durations / precision;
        for i = 1:numel(onsets_adj)
            task_time_course(onsets_adj(i):onsets_adj(i) + duration_adj(i) - 1) = 1;
        end
        % Calculate HRF-convolved task
        hrf = spm_hrf(TR);
        convolved_ttc = spm_Volterra(struct('u', task_time_course, 'name', {{'task'}}), hrf);
        convolved_ttc = convolved_ttc(1:Nvol);
    else
        disp('ERROR: function not working for timing_units = secs');
        %    precision = 0.01;
        %    task_time_course = zeros(Nvol*TR/precision, 1);
        %    onsets_adj = round(onsets, 2)/precision;
        %    duration_adj = round(durations, 2)/precision;
        %    for i = 1:numel(onsets_adj)
        %        task_time_course(onsets_adj(i):onsets_adj(i)+duration_adj(i)-1) = 1;
        %    end
        %    % Calculate HRF-convolved task
        %    hrf = spm_hrf(TR);
        %    convolved_ttc = spm_Volterra(struct('u', task_time_course, 'name', {{'task'}}), hrf);
        %    convolved_ttc = convolved_ttc(1:(Nvol*TR/precision));
    end

    % Condition-encoded task-baseline vector
    % vectEncCond = task_design + 1; % 1 = index in baseline block; 2 = index in task block
    %% Get indices for task and baseline blocks
    % task_blocks = bwlabel(task_design);
    % baseline_blocks = bwlabel(~task_design);
    % baseline_design = (~task_design);
    % I_baseline = [];
    % for n = 1:length(baseline_onsets)
    %    I_block = find(baseline_blocks == n);
    %    if n > 1
    %        I_block = I_block((1+TR_skip):end);
    %    end
    %    I_baseline = [I_baseline I_block];
    % end
