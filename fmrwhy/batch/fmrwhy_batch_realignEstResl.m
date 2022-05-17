function fmrwhy_batch_realignEstResl(functional_fn, template_fn, saveAs_fn, varargin)

    % -------------
    % Parse inputs
    % -------------
    validChar = @(x) ischar(x);
    validType = @(x) any(validatestring(x, filetypes));

    eOptions = {'quality', 'sep', 'fwhm', 'rtm', 'einterp', 'ewrap', 'weight'};
    eDefaults = {0.9, 4, 5, 0, 2, [0 0 0], ''};
    rOptions = {'which', 'rinterp', 'rwrap', 'mask', 'prefix'};
    rDefaults = {[1 0], 4, [0 0 0], 1, 'r'};

    p = inputParser;
    addRequired(p, 'functional_fn');
    addRequired(p, 'template_fn');
    addRequired(p, 'saveAs_fn');
    for i = 1:numel(eOptions)
        addParameter(p, eOptions{i}, eDefaults{i});
    end
    for i = 1:numel(rOptions)
        addParameter(p, rOptions{i}, rDefaults{i});
    end
    parse(p, functional_fn, template_fn, saveAs_fn, varargin{:});
    params = p.Results;
    functional_fn = params.functional_fn;
    template_fn = params.template_fn;
    saveAs_fn = params.saveAs_fn;

    % ------------------------------------------
    % Set up functional filenames for processing
    % ------------------------------------------

    % Could be a single filename, or cell array of filenames
    if iscell(functional_fn)
        N_runs = numel(functional_fn);
        func_files = functional_fn;
        if ~iscell(saveAs_fn) || numel(saveAs_fn) ~= N_runs
            % Add error message TODO
            disp('Error: saveAs_fn should have the same data type and number of indices as functional_fn');
        end

    elseif isstring(functional_fn) || ischar(functional_fn)
        fn_charstring = true;
        N_runs = 1;
        func_files = {functional_fn};
        saveAs_fn = {saveAs_fn}; % TODO first test if the save fn is also a single string/chararray
    else
        % Add error message TODO
        disp('Error: functional file input is not cell array, nor char/string');
    end

    data = cell(1, N_runs);
    temp_functional_fn = cell(1, N_runs);
    for r = 1:N_runs
        [d, f, e] = fileparts(func_files{r});
        temp_functional_fn{r} = fullfile(d, ['temp_' f e]);
        copyfile(func_files{r}, temp_functional_fn{r});
        func_spm = spm_vol(temp_functional_fn{r});
        Nt = numel(func_spm);
        % Filenames for which to estimate 3D realignment parameters
        fns = {};
        if template_fn == 0
            for i = 1:Nt
                fns{i} = [temp_functional_fn{r} ',' num2str(i)];
            end
        else
            fns{1} = [template_fn ',1'];
            for i = 1:Nt
                fns{i + 1} = [temp_functional_fn{r} ',' num2str(i)];
            end
        end
        data{r} = fns';
    end

    % -------------------
    % Set up SPM12 job
    % -------------------

    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    % Data
    realign_estimate_reslice = struct;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.data = data;
    % Eoptions
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = params.quality;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = params.sep;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = params.fwhm;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = params.rtm; % Default = register to first
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = params.einterp;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = params.ewrap;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = params.weight;
    % Roptions
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = params.which; % Refault = images [2..n]
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = params.rinterp;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = params.rwrap;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = params.mask;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = params.prefix;

    % -------------------
    % Run SPM12 job
    % -------------------
    spm_jobman('run', realign_estimate_reslice.matlabbatch);

    % -------------
    % Format output
    % -------------
    output = struct;
    for r = 1:N_runs
        [d, f, e] = fileparts(temp_functional_fn{r});
        rtemp_functional_fn = fullfile(d, ['r' f e]);
        [status, msg, msgID] = movefile(rtemp_functional_fn, saveAs_fn{r});
        delete(temp_functional_fn{r});
    end
