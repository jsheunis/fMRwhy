function confounds_tsv = fmrwhy_preproc_generateMultRegr(bids_dir, sub, ses, task, run, echo)
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

% TODO: %order = 6; % 6 / 12 / 24 ==> MP, MP + MPdiff, MP + MPdiff + MP??2 + MPdiff??2


% ---


% Load/create required defaults
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
template_task = 'motor'; % changed for fingertapping experiment. TODO: change back. and update functioning.
template_run = '1';
template_echo = '2';
stre_txt = ['sub-' sub '_task-' task '_run-' run '_echo-' echo];
str_txt = ['sub-' sub '_task-' task '_run-' run];

% Directory and content setup
deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
qc_dir = fullfile(deriv_dir, 'fmrwhy-qc');
sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);
sub_dir_qc = fullfile(qc_dir, ['sub-' sub]);
func_dir_qc = fullfile(sub_dir_qc, 'func');
if ~exist(func_dir_qc, 'dir')
    mkdir(func_dir_qc)
end

% Defaults to include/exclude
include_fd = 1;
include_tissue = 1;
include_physio = 0;


% Grab functional timeseries filename,
functional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_bold.nii']);

% Grab realigned (and slice time corrected?) functional timeseries filename,
rfunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rpreproc_bold.nii']);
rafunctional_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_bold.nii']);

% Grab template filename
template_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' template_task '_run-' template_run '_space-individual_bold.nii']);


% -------
% STEP 1: Generate realignment, framewise displacement and censoring regressors
% -------
% First realignment parameters
% Check if this has already been done by seeing if the tsv file with head movement parameters exist
motion_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_motion.tsv']);
[d, f, e] = fileparts(motion_fn);
if ~exist(motion_fn, 'file')
    % If it does not exist, estimate MPs
    disp(['Estimating 3D realignment parameters for: ' stre_txt]);
    realign_measures = fmrwhy_batch_realignEst(functional_fn, template_fn);
    temp_txt_fn = fullfile(d, [f '.txt']);
    col_names = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
    MP = realign_measures.MP;
    data_table = array2table(MP,'VariableNames', col_names);
    writetable(data_table, temp_txt_fn, 'Delimiter','\t');
    [status, msg, msgID] = movefile(temp_txt_fn, motion_fn);
    disp('Complete!')
    disp('---')
else
    disp(['3D realignment parameters already estimated: ' motion_fn])
    disp('---')
    MP = struct2array(tdfread(motion_fn));
end
if include_fd
    % Then framewise displacement and censoring
    framewise_displacement_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_fd.tsv']);
    [d, f, e] = fileparts(framewise_displacement_fn);
    if ~exist(framewise_displacement_fn, 'file')
        % If it does not exist, estimate MPs
        disp(['Preparing FD measures for: ' stre_txt]);
        r = 50; % mm
        FD_threshold = 0.5; % mm
        FD_measures = fmrwhy_qc_calculateFD(MP, r, FD_threshold); %FD_measures.FD_outliers_regr
        FD = FD_measures.FD;
        FD_outliers_regr = FD_measures.FD_outliers_regr;
        temp_txt_fn = fullfile(d, [f '.txt']);
        col_names = {'framewise_displacement','framewise_displacement_censor'};
        data_table = array2table([FD FD_outliers_regr],'VariableNames', col_names);
        writetable(data_table, temp_txt_fn, 'Delimiter','\t');
        [status, msg, msgID] = movefile(temp_txt_fn, framewise_displacement_fn);
        disp('Complete!')
        disp('---')
    else
        disp(['FD measures already calculated: ' framewise_displacement_fn])
        disp('---')
        FD_measures = struct2array(tdfread(framewise_displacement_fn));
        FD = FD_measures(:,1);
        FD_outliers_regr = FD_measures(:,2);
    end
end

% -------
% STEP 2: Nuisance signals from tissue compartment and whole brain masks
% TODO: figure out from which preprocessed timeseries this has to be calculated. For now using rafunctional.
% -------
if include_tissue
    tissue_regr_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_tissue.tsv']);
    [d, f, e] = fileparts(tissue_regr_fn);
    if ~exist(tissue_regr_fn, 'file')
        disp(['Calculating tissue compartment signal averages: ' tissue_regr_fn])
        disp('---')
        % Get anatomical tissue compartment masks in individual functional space
        masks = fmrwhy_util_loadMasks(bids_dir, sub);
        % Get timeseries data
        rafunctional_4D = spm_read_vols(spm_vol(rafunctional_fn));
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
        [status, msg, msgID] = movefile(temp_txt_fn, tissue_regr_fn);
    else
        disp(['Tissue compartment signal averages already calculated: ' tissue_regr_fn])
        disp('---')
        tissue_regressors = struct2array(tdfread(tissue_regr_fn));
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
    log_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_physio.tsv.gz']);
    % Create struct with options for PhysIO, TODO: ideally load from workflow settings and/or defaults
    options = struct;
    options.save_dir = fullfile(sub_dir_qc, 'func');;
    options.vendor = 'BIDS';
    options.cardiac_fn = log_fn;
    options.respiration_fn = log_fn;
    options.sampling_interval = 0.002; % 500 Hz ==> Philips wired acquisition
    options.align_scan = 'last';
    options.Nslices = 34;
    options.TR = 2; % in seconds
    options.Ndummies = 5; % include, even if these are not included in the fMRI timeseries data exported from the scanner
    options.Nscans = 210;
    options.onset_slice = 1;
    options.cardiac_modality = 'PPU';
    options.output_multiple_regressors_fn = 'PhysIO_multiple_regressors.txt'; % text file name
    % Run PhysIO if required
    physio_regr_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_physio.tsv']);
    [d, f, e] = fileparts(physio_regr_fn);
    if ~exist(physio_regr_fn, 'file')
        disp('Physio regressors not calculated yet. Calculating now.')
        phys_data = fmrwhy_batch_PhysIO(options);
        temp_txt_fn = fullfile(d, [f '.txt']);
        col_names = {'retroicor_c1','retroicor_c2','retroicor_c3','retroicor_c4','retroicor_c5','retroicor_c6','retroicor_r1','retroicor_r2','retroicor_r3','retroicor_r4','retroicor_r5','retroicor_r6','retroicor_r7','retroicor_r8','retroicor_cxr1','retroicor_cxr2','retroicor_cxr3','retroicor_cxr4','hrv','rvt'};
        data = load(phys_data.multiple_regressors_fn);
        data_table = array2table(data,'VariableNames', col_names);
        writetable(data_table, temp_txt_fn, 'Delimiter','\t');
        [status, msg, msgID] = movefile(temp_txt_fn, physio_regr_fn);
    else
        disp('Physio regressors previously calculated. Loading now.')
        phys_data.physio_regr = struct2array(tdfread(physio_regr_fn));
    end
end


% -------
% STEP 4: Add all regressors to one array and save to tsv
% TODO: this does not replace the txt file, i.e. duplicate text file
% remains. Fix this.
% -------
confounds_fn = fullfile(sub_dir_preproc, 'func', ['sub-' sub '_task-' task '_run-' run '_desc-confounds_regressors.tsv']);
[d, f, e] = fileparts(confounds_fn);
if ~exist(confounds_fn, 'file')
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




