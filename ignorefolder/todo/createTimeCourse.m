function task_time_course = createTimeCourse(Nvol, TR, onsets, durations)

precision = 0.001;
task_time_course = zeros(Nvol*TR/precision, 1);

onsets_adj = onsets/precision;
duration_adj = durations/precision;

for i = 1:numel(onsets_adj)
    task_time_course(onsets_adj(i):onsets_adj(i)+duration_adj(i)-1) = 1;
end
