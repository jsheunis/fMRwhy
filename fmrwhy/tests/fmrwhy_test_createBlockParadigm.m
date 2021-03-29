% A custom workflow that does ...


%--------------------------------------------------------------------------

% Load/create required parameters
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Loop through subjects, sessions, tasks, runs, etc
sub = '001';
ses = '';
task = 'emotion';
run = '1';
echo = '2';

events_fn = fullfile(bids_dir, 'derivatives', 'fmrwhy-preproc', ['sub-' sub], 'func',  ['sub-' sub '_task-' task '_run-' run '_events.tsv']);
cond_names = {'Faces', 'Shapes'};

[cond, trials] = fmrwhy_util_1stlevelBIDStoConditions(events_fn, cond_names);

TR = 2;
timing_units = 'secs';
onsets = cond(1).onset;
durations = cond(1).duration;
precision = 0.01;
Nvol = 210;



[task_time_course, convolved_ttc, hrf] = fmrwhy_util_createBlockParadigm(Nvol, TR, onsets, durations, precision, timing_units);


figure;
ax = subplot(2,1,1);
im_task = imagesc(ax, task_time_course'); colormap(gray);
xlim(ax,[0 Nvol*TR/precision]);
ylim(ax,[0.5 1.55]);

ax2 = subplot(2,1,2);
plot(convolved_ttc);
xlim(ax2,[0 Nvol*TR/precision]);