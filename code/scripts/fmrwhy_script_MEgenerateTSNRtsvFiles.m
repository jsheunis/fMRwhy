% -------
% STEP 1 -- Load defaults, filenames and parameters
% -------
% Load fMRwhy defaults
options = fmrwhy_defaults;

% Main input: BIDS root folder
bids_dir = '/Volumes/Stephan_WD/NEUFEPME_data_BIDS';

% Setup fmrwhy BIDS-derivatuve directories on workflow level
options = fmrwhy_defaults_setupDerivDirs(bids_dir, options);

% Grab parameters from workflow settings file
options = fmrwhy_settings_preprocQC(bids_dir, options);

options.dash_me_dir = '/Users/jheunis/Documents/Websites/rt-me-fmri-dash/bids/derivatives/fmrwhy-multiecho';
options.dash_deriv_dir = '/Users/jheunis/Documents/Websites/rt-me-fmri-dash/bids/derivatives';

% Set subject, sessions
subs = {'002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
for s = 1:numel(subs)
    sub = subs{s};
    ses = '';

    % Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
    options = fmrwhy_defaults_setupSubDirs(bids_dir, sub, options);

    % Update workflow params with subject anatomical derivative filenames
    options = fmrwhy_defaults_subAnat(bids_dir, sub, options);

    % load mask
    masks = fmrwhy_util_loadOrientMasks(bids_dir, sub);
    mask_fn = masks.brain_mask_fn;

    options.me_dir = fullfile(options.deriv_dir, 'fmrwhy-multiecho');
    options.sub_dir_me = fullfile(options.me_dir, ['sub-' sub]);
    options.func_dir_me = fullfile(options.sub_dir_me, 'func');
    dash_sub_dir = fullfile(options.dash_me_dir, ['sub-' sub]);
    if ~exist(dash_sub_dir, 'dir')
        mkdir(dash_sub_dir)
    end

    % Loop through sessions, tasks, runs, etc
    tasks = {'rest', 'motor', 'emotion'};
    runs = {'1', '2'};
    echo = '2';


    for t = 1:numel(tasks)

        task = tasks{t};

        for r = 1:numel(runs)
            run = runs{r};

            disp('------------')
            disp('------------')
            disp(['Task: ' task ';  Run: ' run])
            disp('------------')
            disp('------------')

            % Filenames
            options = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, options);

            % Calculate tSNR for each timeseries
            rafunctional_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_echo-2_desc-rapreproc_tsnr.nii']);
            combined_t2s_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEt2star_tsnr.nii']);
            combined_tsnr_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEtsnr_tsnr.nii']);
            combined_te_fn = fullfile(options.func_dir_me, ['sub-' sub '_task-' task '_run-' run '_desc-combinedMEte_tsnr.nii']);
            tsnr_fns = {rafunctional_fn, combined_t2s_fn, combined_tsnr_fn, combined_te_fn};
            tsnr_output_fns = {};
            for i = 1:numel(tsnr_fns)

                if strcmp(task, 'rest') == 1 && strcmp(run, '1') == 1 && i > 1
                    disp('------------')
                    disp(['... Skipping combined echoes for task: ' task ';  Run: ' run ' ...'])
                    disp('------------')
                    continue;
                end

                [p_tsnr, frm, rg, dim] = fmrwhy_util_readOrientNifti(tsnr_fns{i});
                tsnr_img = p_tsnr.nii.img(:);
    %            fmrwhy_util_calculateTSNR(main_fns{i}, mask_fn, tsnr_fns{i}, template_fn);
                for j = 1:4
                    vals = tsnr_img(masks.([masks.field_names{j} '_mask_I']));
                    tsnr_output_fns{i,j} = strrep(tsnr_fns{i}, '_tsnr.nii', ['_' masks.field_names{j} 'tsnr.tsv']);

                    temp_txt_fn = strrep(tsnr_output_fns{i,j}, '.tsv', '_temp.txt');

                    data_table = array2table(vals,'VariableNames', {'tsnr'});
                    writetable(data_table, temp_txt_fn, 'Delimiter','\t');
                    [status, msg, msgID] = movefile(temp_txt_fn, tsnr_output_fns{i,j});

                    [d,f,e] = fileparts(tsnr_output_fns{i,j});
                    new_fn = fullfile(dash_sub_dir, [f e]);
                    copyfile(tsnr_output_fns{i,j}, new_fn)
                end
            end
        end
    end
end