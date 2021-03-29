% -------
% STEP 1 -- Load defaults, filenames and parameters
% -------
% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

% Set subject, sessions
sub = '001';
ses = '';

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

% Update workflow params with subject anatomical derivative filenames
options = fmrwhy_defaults_subAnat(bids_dir, sub, options);




for e = 1:options.Ne
    disp('---')
    disp(['Echo ' num2str(e)])
    disp('---')
    % Update workflow params with subject functional derivative filenames
    options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, num2str(e), options);
    tsnr_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' num2str(e) '_desc-rapreproc_tsnr.nii']);

    disp(['Temporal SNR file does not exist yet. Calculating now...']);
    tsnr_output = fmrwhy_util_calculateTSNR(options.rafunctional_fn, 0, tsnr_fn, options.template_fn)
end



% Loop through sessions, tasks, runs, etc
tasks = {'rest'};
runs = {'2'};
echo = '2';




for t = 1:numel(tasks)

    task = tasks{t};

    for r = 1:numel(runs)
        run = runs{r};

        if strcmp(task, 'rest') == 1 && strcmp(run, '1') == 1
            disp('------------')
            disp(['... Skipping Task: ' task ';  Run: ' run ' ...'])
            disp('------------')
            continue;
        end

        disp('------------')
        disp('------------')
        disp(['Task: ' task ';  Run: ' run])
        disp('------------')
        disp('------------')

        % Filenames
        options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

        t2star_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-1_desc-MEparams_t2star.nii']);
        s0_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-1_desc-MEparams_s0.nii']);
        combined_t2s_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEt2star_bold.nii']);
        combined_tsnr_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEtsnr_bold.nii']);
        combined_te_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEte_bold.nii']);
        combined_fns = {combined_t2s_fn, combined_tsnr_fn, combined_te_fn};

        % Grab and construct parameters (data and weights) for multi-echo combination
        disp('Concatenating functional data')
        TE = options.TE;
        template_spm = spm_vol(options.template_fn);
        template_dim = template_spm.dim;
        Nt = options.Nscans;
        sz = [template_dim Nt numel(TE)];
        func_data = zeros(sz);
        for e = 1:numel(TE)
            echo = num2str(e);
            rafunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-' echo '_desc-rapreproc_bold.nii'])
            func_data(:,:,:,:,e) = spm_read_vols(spm_vol(rafunctional_fn));
%            nii = nii_tool('load', rafunctional_fn);
%            func_data(:,:,:,:,e) = double(nii.img);
        end

        disp('Loading weight images')
        t2star_img = spm_read_vols(spm_vol(t2star_fn));
%        nii = nii_tool('load', t2star_fn);
%        t2star_img = double(nii.img);
        tsnr_data = zeros([template_dim numel(TE)]);
        for e = 1:numel(TE)
            echo = num2str(e);
            tsnr_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_echo-' echo '_desc-rapreproc_tsnr.nii']);
            tsnr_data(:,:,:,e) = spm_read_vols(spm_vol(tsnr_fn));
%            nii = nii_tool('load', tsnr_fn);
%            tsnr_data(:,:,:,e) = double(nii.img);
        end
        % Mult-echo combination
        disp('Combining multiple echoes with different combination methods')
        combined_dataAll_t2s = fmrwhy_me_combineEchoes(func_data, TE, 0, 1, t2star_img);
        combined_dataAll_tsnr = fmrwhy_me_combineEchoes(func_data, TE, 0, 2, tsnr_data);
        combined_dataAll_TE = fmrwhy_me_combineEchoes(func_data, TE, 0, 3, TE);

        % Save nifti images
        % TODO: replace this with dicm2nii equivalent (ask on github repo)
        rafunctional_fn = fullfile(options.func_dir_preproc, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_bold.nii']);
        new_spm_t2s = spm_vol(rafunctional_fn);
        new_spm_tsnr = new_spm_t2s;
        new_spm_te = new_spm_t2s;

        for i = 1:numel(new_spm_t2s)
            new_spm_t2s(i).fname = combined_t2s_fn;
            new_spm_t2s(i).private.dat.fname = combined_t2s_fn;
            spm_write_vol(new_spm_t2s(i), combined_dataAll_t2s(:,:,:,i));

            new_spm_tsnr(i).fname = combined_tsnr_fn;
            new_spm_tsnr(i).private.dat.fname = combined_tsnr_fn;
            spm_write_vol(new_spm_tsnr(i), combined_dataAll_tsnr(:,:,:,i));

            new_spm_te(i).fname = combined_te_fn;
            new_spm_te(i).private.dat.fname = combined_te_fn;
            spm_write_vol(new_spm_te(i), combined_dataAll_TE(:,:,:,i));
        end
    end
end