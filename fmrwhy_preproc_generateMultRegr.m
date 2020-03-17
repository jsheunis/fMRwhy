function confounds_tsv = fmrwhy_preproc_generateMultRegr(bids_dir, sub, ses, task, run, echo, opts)
% For now, this function assumes that all necessary preproc has been done and all necessary files exist
% Generates multiple regressors for quality control and

% trans_x trans_y trans_z rot_x rot_y rot_z framewise_displacement framewise_displacement_censor global_signal white_matter csf
% std_dvars dvars trans_x_derivative1 trans_x_power2 ....
% retroicor_c1 retroicor_c2 retroicor_c3 retroicor_c4 retroicor_c5 retroicor_c6
% retroicor_r1 retroicor_r2 retroicor_r3 retroicor_r4 retroicor_r5 retroicor_r6 retroicor_r7 retroicor_r8
% hrv rvt

%RETROICOR cardiac regressors [2 x nOrderCardiac] ==> 2 x 3 = 6
%RETROICOR respiratory regressors [2 x nOrderRespiratory] ==> 2 x 4 = 8
%RETROICOR cardXResp interaction regressors [4 x nOrderCardiacXRespiratory] ==> 4 x 1 = 4
%HRV [nDelaysHRV]
%RVT [nDelaysRVT]

% trans_x_derivative1 trans_y_derivative1 trans_z_derivative1 rot_x_derivative1 rot_y_derivative1 rot_z_derivative1
% trans_x_power2 trans_y_power2 trans_z_power2 rot_x_power2 rot_y_power2 rot_z_power2
% trans_x_derivative1_power2 trans_y_derivative1_power2 trans_z_derivative1_power2 rot_x_derivative1_power2 rot_y_derivative1_power2 rot_z_derivative1_power2

% TODO: %order = 6; % 6 / 12 / 24 ==> MP, MP + MPdiff, MP + MPdiff + MP??2 + MPdiff??2


% ---


% Setup fmrwhy bids directories on workflow level
fmrwhy_defaults_setupDerivDirs(bids_dir);

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
fmrwhy_defaults_setupSubDirs(bids_dir, sub);

% Update workflow params with subject anatomical derivative filenames
opts = fmrwhy_defaults_subAnat(bids_dir, sub, opts);

% Update workflow params with subject functional derivative filenames
opts = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, opts);

% Filler text
stre_txt = ['sub-' sub '_task-' task '_run-' run '_echo-' echo];
str_txt = ['sub-' sub '_task-' task '_run-' run];

% Default steps to include/exclude
include_fd = opts.confounds.include_fd;
include_tissue = opts.confounds.include_tissue;
include_physio = opts.confounds.include_physio;

% -------
% STEP 1: Generate realignment, framewise displacement and censoring regressors
% -------
% First realignment parameters
% Check if this has already been done by seeing if the tsv file with head movement parameters exist
[d, f, e] = fileparts(opts.motion_fn);
if ~exist(opts.motion_fn, 'file')
    % If it does not exist, estimate MPs
    disp(['Estimating 3D realignment parameters for: ' stre_txt]);
    realign_measures = fmrwhy_batch_realignEst(opts.functional_fn, opts.template_fn);
    temp_txt_fn = fullfile(d, [f '.txt']);
    col_names = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
    MP = realign_measures.MP;
    data_table = array2table(MP,'VariableNames', col_names);
    writetable(data_table, temp_txt_fn, 'Delimiter','\t');
    [status, msg, msgID] = movefile(temp_txt_fn, opts.motion_fn);
    disp('Complete!')
    disp('---')
else
    disp(['3D realignment parameters already estimated: ' opts.motion_fn])
    disp('---')
    MP = struct2array(tdfread(opts.motion_fn));
end

if include_fd
    % Then framewise d isplacement and censoring
    [d, f, e] = fileparts(opts.framewise_displacement_fn);
    if ~exist(opts.framewise_displacement_fn, 'file')
        % If it does not exist, estimate MPs
        disp(['Preparing FD measures for: ' stre_txt]);
        FD_measures = fmrwhy_qc_calculateFD(MP, opts.r, opts.FD_threshold); %FD_measures.FD_outliers_regr
        FD = FD_measures.FD;
        FD_outliers_regr = FD_measures.FD_outliers_regr;
        temp_txt_fn = fullfile(d, [f '.txt']);
        col_names = {'framewise_displacement','framewise_displacement_censor'};
        data_table = array2table([FD FD_outliers_regr],'VariableNames', col_names);
        writetable(data_table, temp_txt_fn, 'Delimiter','\t');
        [status, msg, msgID] = movefile(temp_txt_fn, opts.framewise_displacement_fn);
        disp('Complete!')
        disp('---')
    else
        disp(['FD measures already calculated: ' opts.framewise_displacement_fn])
        disp('---')
        FD_measures = struct2array(tdfread(opts.framewise_displacement_fn));
        FD = FD_measures(:,1);
        FD_outliers_regr = FD_measures(:,2);
    end
end

% -------
% STEP 2: Nuisance signals from tissue compartment and whole brain masks
% TODO: figure out from which preprocessed timeseries this has to be calculated. For now using rafunctional.
% -------
if include_tissue
    [d, f, e] = fileparts(opts.tissue_regr_fn);
    if ~exist(opts.tissue_regr_fn, 'file')
        disp(['Calculating tissue compartment signal averages: ' opts.tissue_regr_fn])
        disp('---')
        % Get anatomical tissue compartment masks in individual functional space
        masks = fmrwhy_util_loadMasks(bids_dir, sub);
        % Get timeseries data
        rafunctional_4D = spm_read_vols(spm_vol(opts.rafunctional_fn));
        [Ni, Nj, Nk, Nt] = size(rafunctional_4D);
        rafunctional_2D = reshape(rafunctional_4D, Nj*Nj*Nk, Nt); % [Nj*Nj*Nk, Nt]
        % Get spatial mean per timepoint per tissue compartment
        GM_timeseries = (mean(rafunctional_2D(masks.GM_mask_I, :)))';
        WM_timeseries = (mean(rafunctional_2D(masks.WM_mask_I, :)))';
        CSF_timeseries = (mean(rafunctional_2D(masks.CSF_mask_I, :)))';
        brain_timeseries = (mean(rafunctional_2D(masks.brain_mask_I, :)))';
        % Save to tsv file
        temp_txt_fn = fullfile(d, [f '.txt']);
        col_names = {'grey_matter', 'white_matter', 'csf', 'global_signal'};
        data = [GM_timeseries WM_timeseries CSF_timeseries brain_timeseries];
        data_table = array2table(data,'VariableNames', col_names);
        writetable(data_table, temp_txt_fn, 'Delimiter','\t');
        [status, msg, msgID] = movefile(temp_txt_fn, opts.tissue_regr_fn);
    else
        disp(['Tissue compartment signal averages already calculated: ' opts.tissue_regr_fn])
        disp('---')
        tissue_regressors = struct2array(tdfread(opts.tissue_regr_fn));
        GM_timeseries = tissue_regressors(:,1);
        WM_timeseries = tissue_regressors(:,2);
        CSF_timeseries = tissue_regressors(:,3);
        brain_timeseries = tissue_regressors(:,4);
    end
end

% -------
% STEP 3: Physiological regressors (Retroicor, HRV, RVT)
% -------
if include_physio
    % Grab physio log-file, e.g. sub-001_task-rest_run-1_physio.tsv.gz
    log_fn = fullfile(opts.sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_physio.tsv.gz']);

    % Create struct with options for PhysIO
    physio_options = opts.physio.options;
    physio_options.save_dir = fullfile(sub_dir_qc, 'func');;
    physio_options.cardiac_fn = log_fn;
    physio_options.respiration_fn = log_fn;

    % Run PhysIO if required
    [d, f, e] = fileparts(opts.physio_regr_fn);
    if ~exist(opts.physio_regr_fn, 'file')
        disp('Physio regressors not calculated yet. Calculating now.')
        phys_data = fmrwhy_batch_PhysIO(physio_options);
        temp_txt_fn = fullfile(d, [f '.txt']);
        col_names = {'retroicor_c1','retroicor_c2','retroicor_c3','retroicor_c4','retroicor_c5','retroicor_c6','retroicor_r1','retroicor_r2','retroicor_r3','retroicor_r4','retroicor_r5','retroicor_r6','retroicor_r7','retroicor_r8','retroicor_cxr1','retroicor_cxr2','retroicor_cxr3','retroicor_cxr4','hrv','rvt'};
        data = load(phys_data.multiple_regressors_fn);
        data_table = array2table(data,'VariableNames', col_names);
        writetable(data_table, temp_txt_fn, 'Delimiter','\t');
        [status, msg, msgID] = movefile(temp_txt_fn, opts.physio_regr_fn);
    else
        disp('Physio regressors previously calculated. Loading now.')
        phys_data.physio_regr = struct2array(tdfread(opts.physio_regr_fn));
    end
end


% -------
% STEP 4: Add all regressors to one array and save to tsv
% TODO: this does not replace the txt file, i.e. duplicate text file
% remains. Fix this.
% -------
[d, f, e] = fileparts(opts.confounds_fn);
if ~exist(opts.confounds_fn, 'file')
    disp('Main confound regressor file not created yet. Creating now.')
    confounds_txt = fullfile(d, [f '.txt']);
    col_names = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
    data = MP;
    if include_fd
        col_names = [col_names, {'framewise_displacement','framewise_displacement_censor'}];
        data = [data FD FD_outliers_regr];
    end
    if include_tissue
        col_names = [col_names, {'white_matter', 'csf', 'global_signal'}];
        data = [data WM_timeseries CSF_timeseries brain_timeseries];
    end
    if include_physio
        col_names = [col_names, {'retroicor_c1','retroicor_c2','retroicor_c3','retroicor_c4','retroicor_c5','retroicor_c6','retroicor_r1','retroicor_r2','retroicor_r3','retroicor_r4','retroicor_r5','retroicor_r6','retroicor_r7','retroicor_r8','retroicor_cxr1','retroicor_cxr2','retroicor_cxr3','retroicor_cxr4','hrv','rvt'}];
        data = [data phys_data.physio_regr];
    end
    dlmwrite(confounds_txt, data, 'delimiter', '\t', 'precision', '%1.7e')
    confounds_tsv = fmrwhy_util_saveAsTSV(confounds_txt, col_names);
else
    disp('Main confound regressor file already exists.')
end




