% Load/create required parameters
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, []);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Loop through subjects, sessions, tasks, runs, etc
sub = '001';

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

% Loop through sessions, tasks, runs, etc
ses = '';
task = 'rest';
run = '1';
echo = '2';

%% Update workflow params with subject functional derivative filenames
options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

include_physio = true;
% -------
% STEP 3: Physiological regressors (Retroicor, HRV, RVT)
% -------
if include_physio
    % Grab physio log-file, e.g. sub-001_task-rest_run-1_physio.tsv.gz
    log_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_physio.tsv.gz']);

    % Create struct with options for PhysIO
    physio_options = options.physio.options;
    physio_options.save_dir = fullfile(options.sub_dir_qc, 'func', ['PhysIO_task-' task '_run-' run]);
    physio_options.cardiac_fn = log_fn;
    physio_options.respiration_fn = log_fn;
    physio_options.level = 0; % verbose.level = 0 ==> do not generate figure outputs during batch process

    % Run PhysIO if required
    [d, f, e] = fileparts(options.physio_regr_fn);
    if ~exist(options.physio_regr_fn, 'file')
        disp('Physio regressors not calculated yet. Calculating now.');
        % Run batch process without generating figures (batch calls `tapas_physio_main_create_regressors`)
        phys_data = fmrwhy_batch_PhysIO(physio_options);
        temp_txt_fn = fullfile(d, [f '.txt']);
        col_names = {'retroicor_c1', 'retroicor_c2', 'retroicor_c3', 'retroicor_c4', 'retroicor_c5', 'retroicor_c6', 'retroicor_r1', 'retroicor_r2', 'retroicor_r3', 'retroicor_r4', 'retroicor_r5', 'retroicor_r6', 'retroicor_r7', 'retroicor_r8', 'retroicor_cxr1', 'retroicor_cxr2', 'retroicor_cxr3', 'retroicor_cxr4', 'hrv', 'rvt'};
        data = load(phys_data.multiple_regressors_fn);
        data_table = array2table(data, 'VariableNames', col_names);
        writetable(data_table, temp_txt_fn, 'Delimiter', '\t');
        [status, msg, msgID] = movefile(temp_txt_fn, options.physio_regr_fn);
        % Generate figures (by calling `tapas_physio_review` with updated physio options)
        phys_file = fullfile(physio_options.save_dir, 'physio.mat');
        physmat = load(phys_file);
        physmat.physio.verbose.level = 2;
        physmat.physio.verbose.fig_output_file = ['sub-' sub '_task-' task '_run-' run '_physioQC.jpg'];
        physmat.physio.verbose.show_figs = false;
        physmat.physio.verbose.save_figs = true;
        physmat.physio.verbose.close_figs = true;
        % Run tapas_physio_review
        phys_data_figures = tapas_physio_review(physmat.physio);

    else
        disp('Physio regressors previously calculated. Loading now.');
        phys_data.physio_regr = struct2array(tdfread(options.physio_regr_fn));
    end
end

% TAPAS
%% Run PhysIO if required
%    [d, f, e] = fileparts(options.physio_regr_fn);
%    if ~exist(options.physio_regr_fn, 'file')
%        disp('Physio regressors not calculated yet. Calculating now.')
%        % Run batch process without generating figures (batch calls `tapas_physio_main_create_regressors`)
%        phys_data = fmrwhy_tapas_PhysIO(physio_options);
%        temp_txt_fn = fullfile(d, [f '.txt']);
%        col_names = {'retroicor_c1','retroicor_c2','retroicor_c3','retroicor_c4','retroicor_c5','retroicor_c6','retroicor_r1','retroicor_r2','retroicor_r3','retroicor_r4','retroicor_r5','retroicor_r6','retroicor_r7','retroicor_r8','retroicor_cxr1','retroicor_cxr2','retroicor_cxr3','retroicor_cxr4','hrv','rvt'};
%        data = phys_data.R;
%        data_table = array2table(data,'VariableNames', col_names);
%        writetable(data_table, temp_txt_fn, 'Delimiter','\t');
%        [status, msg, msgID] = movefile(temp_txt_fn, options.physio_regr_fn);
%        % Generate figures (by calling `tapas_physio_review` with updated physio options)
%        physio = phys_data.physio_mat;
%        physio.verbose.level = 2;
%        physio.verbose.fig_output_file = ['sub-' sub '_task-' task '_run-' run '_physioQC.jpg'];
%        physio.verbose.show_figs = false;
%        physio.verbose.save_figs = true;
%        physio.verbose.close_figs = true;
%        % Run tapas_physio_review
%        phys_data_figures = tapas_physio_review(physio);
%
%    else
%        disp('Physio regressors previously calculated. Loading now.')
%        phys_data.physio_regr = struct2array(tdfread(options.physio_regr_fn));
%    end

% BATCH
%% Run PhysIO if required
%    [d, f, e] = fileparts(options.physio_regr_fn);
%    if ~exist(options.physio_regr_fn, 'file')
%        disp('Physio regressors not calculated yet. Calculating now.')
%        % Run batch process without generating figures (batch calls `tapas_physio_main_create_regressors`)
%        phys_data = fmrwhy_batch_PhysIO(physio_options);
%        temp_txt_fn = fullfile(d, [f '.txt']);
%        col_names = {'retroicor_c1','retroicor_c2','retroicor_c3','retroicor_c4','retroicor_c5','retroicor_c6','retroicor_r1','retroicor_r2','retroicor_r3','retroicor_r4','retroicor_r5','retroicor_r6','retroicor_r7','retroicor_r8','retroicor_cxr1','retroicor_cxr2','retroicor_cxr3','retroicor_cxr4','hrv','rvt'};
%        data = load(phys_data.multiple_regressors_fn);
%        data_table = array2table(data,'VariableNames', col_names);
%        writetable(data_table, temp_txt_fn, 'Delimiter','\t');
%        [status, msg, msgID] = movefile(temp_txt_fn, options.physio_regr_fn);
%        % Generate figures (by calling `tapas_physio_review` with updated physio options)
%        phys_file = fullfile(physio_options.save_dir, 'physio.mat');
%        physmat = load(phys_file);
%        physmat.physio.verbose.level = 2;
%        physmat.physio.verbose.fig_output_file = ['sub-' sub '_task-' task '_run-' run '_physioQC.jpg'];
%        physmat.physio.verbose.show_figs = false;
%        physmat.physio.verbose.save_figs = true;
%        physmat.physio.verbose.close_figs = true;
%        % Run tapas_physio_review
%        phys_data_figures = tapas_physio_review(physmat.physio);
%
%    else
%        disp('Physio regressors previously calculated. Loading now.')
%        phys_data.physio_regr = struct2array(tdfread(options.physio_regr_fn));
%    end
