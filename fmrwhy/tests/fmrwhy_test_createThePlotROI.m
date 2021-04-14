options = fmrwhy_defaults;
% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';
% Loop through subjects, sessions, tasks, runs, etc
sub = '001';
% Loop through sessions, tasks, runs, etc
ses = '';
task = 'emotion';
run = '2';
echo = '2';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% Update workflow params with subject functional derivative filenames
options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

% Ignore the 'rest' task (assume there is no task ROI for this; have to change in future if RSnetworks available to be normalised or something)
if strcmp(task, 'rest') ~= 1
    % Loop through all ROIs for the particular task
    for j = 1:numel(options.roi.(task).orig_fn)
        functional_fn = options.sfunctional_fn;
        desc = options.roi.(task).desc{j};
        saveAs_fn = fullfile(options.func_dir_qc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-' desc '_grayplot.png']);

        task_info.TR = options.firstlevel.(task).sess_params.timing_RT;
        task_info.onsets = options.firstlevel.(task).sess_params.cond_onset;
        task_info.durations = options.firstlevel.(task).sess_params.cond_duration;
        task_info.precision = 1;

        if ~exist(saveAs_fn, 'file')
            trace_info = [];
            fmrwhy_util_thePlotROI(functional_fn, options.brain_mask_fn, options.roi.(task).rroi_fn{j}, task_info, trace_info, saveAs_fn);
        else
            disp(['File already exists: ' saveAs_fn]);
        end
    end
else
    disp('---');
    disp('Not creating ROI timeseries plots for task = rest.');
    disp('---');
end
