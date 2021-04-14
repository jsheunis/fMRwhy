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

% Load multiple confound regressors
confounds_struct = tdfread(options.confounds_fn);
confounds_mat = struct2array(confounds_struct);

physmat_fn = fullfile(options.func_dir_qc, ['PhysIO_task-' task '_run-' run], 'physio.mat');
physio = load(physmat_fn);
card = physio.physio.ons_secs.c;
resp = 2 + physio.physio.ons_secs.r;
% resp = fmrwhy_util_scale(resp,-0.5,0.5);
t_phys = 1:numel(card);

Nt = options.Nscans;

f = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'visible', 'on');

ax2 = subplot(3, 1, 1);
plot(ax2, confounds_struct.trans_x', 'LineWidth', 2);
hold(ax2, 'on');
plot(ax2, confounds_struct.trans_y', 'LineWidth', 2);
plot(ax2, confounds_struct.trans_z', 'LineWidth', 2);
hold(ax2, 'off');
xlim(ax2, [0 Nt]);
ylim(ax2, [-1 1]);
set(ax2, 'Xticklabel', []);
set(ax2, 'Yticklabel', []);
ax2.XAxis.Visible = 'off';
ax2.YAxis.Visible = 'off';

% Ax3 - Physiology
ax3 = subplot(3, 1, 2);
plot(ax3, resp, 'LineWidth', 2);
hold(ax3, 'on');
plot(ax3, card, 'LineWidth', 2);
hold(ax3, 'off');
xlim(ax3, [0 numel(card)]);
ylim(ax3, [-1 4]);
set(ax3, 'Xticklabel', []);
set(ax3, 'Yticklabel', []);
ax3.XAxis.Visible = 'off';
ax3.YAxis.Visible = 'off';
% fmrwhy_util_stretchAx(ax3)
% ax3.Position(1) = axpos(1); ax3.Position(3) = axpos(3);

% Ax4 - Framewise displacement
ax4 = subplot(3, 1, 3);
plot(ax4, confounds_struct.framewise_displacement, 'LineWidth', 2);
xlim(ax4, [0 Nt]);
ylim(ax4, [0 1]);
set(ax4, 'Xticklabel', []);
set(ax4, 'Yticklabel', []);
ax4.XAxis.Visible = 'off';
ax4.YAxis.Visible = 'off';
% fmrwhy_util_stretchAx(ax4)
% ax4.Position(1) = axpos(1); ax4.Position(3) = axpos(3);
