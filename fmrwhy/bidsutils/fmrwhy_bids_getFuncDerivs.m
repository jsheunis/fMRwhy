function options = fmrwhy_bids_getFuncDerivs(bids_dir, sub, task, options, varargin)

    % minimal required arguments:
    %  - positional arguments
    %  - options.preproc_dir
    %  - options.qc_dir

    % -------------
    % Parse inputs
    % -------------
    filetypes = {'func'};
    descriptions = {'Session', 'Acquisition', 'Contrast Enhancing Agent', 'Reconstruction', 'Phase-Encoding Direction', 'Run', 'Echo'};
    % entities = {'ses', 'acq', 'ce', 'rec', 'dir', 'run', 'echo'}; % these entities are required/optional for func bold data specifically (not other types!)
    % looks like bids-matlab is not adhering to bids v1.4.0 and is not able to parse 'ce' and 'dir' entity.
    % for now, remove entities from fMRwhy; have to update bids-matlab in future.
    entities = {'ses', 'acq', 'rec', 'run', 'echo'}; % these entities are required/optional for func bold data specifically (not other types!)
    % {'sub', 'ses', 'acq', 'ce', 'rec', 'fa', 'echo', 'inv', 'run'} % bids-matlab list seems to be pre-v1.4.0
    formats = {'label', 'label', 'label', 'label', 'label', 'index', 'index'};

    validChar = @(x) ischar(x);
    validType = @(x) any(validatestring(x, filetypes));

    p = inputParser;
    addRequired(p, 'bids_dir', validChar);
    addRequired(p, 'sub', validChar);
    addRequired(p, 'task', validChar);
    addRequired(p, 'options');
    for i = 1:numel(entities)
        addParameter(p, entities{i}, '', validChar);
    end
    parse(p, bids_dir, sub, task, options, varargin{:});
    params = p.Results;

    options = params.options;

    % Template filename
    % options.template_fn = fullfile(options.sub_dir_preproc, 'func', ['sub-' sub '_task-' options.template_task '_run-' options.template_run '_space-individual_bold.nii']);

    % Outputs from basicFunc processing

    fields_w_echo = {'motion_fn', 'functional_fn', 'afunctional_fn', 'rfunctional_fn', 'rafunctional_fn', 'sfunctional_fn', 'srfunctional_fn', 'srafunctional_fn'};
    fields_w_echo_desc = {'confounds', '', 'apreproc', 'rpreproc', 'rapreproc', 'spreproc', 'srpreproc', 'srapreproc'};
    fields_w_echo_ext = {'_motion.tsv', '_bold.nii', '_bold.nii', '_bold.nii', '_bold.nii', '_bold.nii', '_bold.nii', '_bold.nii'};

    fields_wo_echo = {'framewise_displacement_fn', 'tissue_regr_fn', 'physio_regr_fn', 'confounds_fn'};
    fields_wo_echo_ext = {'_fd.tsv', '_tissue.tsv', '_physio.tsv', '_regressors.tsv'};

    fields_statsqc = {'mean_fn', 'std_fn', 'tsnr_fn', 'var_fn'};
    fields_statsqc_tsv = {'stats_timeseries_fn', 'stats_summary_fn'};
    fields_statsqc_ext = {'_mean.nii', '_std.nii', '_tsnr.nii', '_var.nii'};

    % Fields with echo included (if multiecho or not)
    [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', params.sub, 'ses', params.ses, 'task', params.task, 'acq', params.acq, 'rec', params.rec, 'run', params.run, 'echo', params.echo);
    options.current_functional_filename = filename;
    for i = 1:numel(fields_w_echo)
        if isempty(fields_w_echo_desc{i})
            options.(fields_w_echo{i}) = fullfile(options.preproc_dir, filepath, [filename fields_w_echo_ext{i}]);
        else
            options.(fields_w_echo{i}) = fullfile(options.preproc_dir, filepath, [filename '_desc-' fields_w_echo_desc{i} fields_w_echo_ext{i}]);
        end

    end

    % Fields with echo excluded (if multiecho or not)
    [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', params.sub, 'ses', params.ses, 'task', params.task, 'acq', params.acq, 'rec', params.rec, 'run', params.run);
    for i = 1:numel(fields_wo_echo)
        options.(fields_wo_echo{i}) = fullfile(options.preproc_dir, filepath, [filename '_desc-confounds' fields_wo_echo_ext{i}]);
    end

    % Stats QC fields with echo excluded (if multiecho or not)
    [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', params.sub, 'ses', params.ses, 'task', params.task, 'acq', params.acq, 'rec', params.rec, 'run', params.run, 'space', 'individual');
    for i = 1:numel(fields_statsqc)
        options.(fields_statsqc{i}) = fullfile(options.qc_dir, filepath, [filename fields_statsqc_ext{i}]);
    end

    for i = 1:numel(fields_statsqc_tsv)
        [filename, filepath] = fmrwhy_bids_constructFilename('func', 'sub', params.sub, 'ses', params.ses, 'task', params.task, 'acq', params.acq, 'rec', params.rec, 'run', params.run, 'desc', fields_statsqc_tsv{i}(1:end - 3), 'ext', '.tsv');
        options.(fields_statsqc_tsv{i}) = fullfile(options.qc_dir, filepath, filename);
    end

    options.basic_func_out_fns = {options.motion_fn, options.afunctional_fn, options.rfunctional_fn, options.rafunctional_fn, options.sfunctional_fn, options.srfunctional_fn, options.srafunctional_fn, options.confounds_fn};
    options.stats_qc_out_fns = {options.mean_fn, options.std_fn, options.tsnr_fn, options.var_fn, options.stats_timeseries_fn, options.stats_summary_fn};
