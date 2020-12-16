function [cond, trials] = fmrwhy_util_1stlevelBIDStoConditions(events_fn, cond_names)


events = tdfread(events_fn);
[r, c] = size(events.trial_type);
trials = zeros(r,1);
for i = 1:numel(cond_names)
    for j = 1:r
        if strcmp(strtrim(events.trial_type(j,:)), cond_names{i})
            trials(j) = i;
        end
    end
    ind = trials == i;
    cond(i).name = cond_names{i};
    cond(i).onset = events.onset(ind);
    cond(i).duration = events.duration(ind);
    cond(i).tmod = 0;
    cond(i).pmod = {''};
    cond(i).orth = 1;
end