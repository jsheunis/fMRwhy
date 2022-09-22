function fmrwhy_batch_VDMApply(functional_fn, vdm_fn, saveAs_fn, varargin)
     

     % -------------
    % Parse inputs
    % -------------
    validChar = @(x) ischar(x);
    validType = @(x) any(validatestring(x, filetypes));

    p = inputParser;
    addRequired(p, 'functional_fn');
    addRequired(p, 'vdm_fn');
    addRequired(p, 'saveAs_fn');
    parse(p, functional_fn, vdm_fn, saveAs_fn, varargin{:});
    params = p.Results;
    functional_fn = params.functional_fn;
    vdm_fn = params.vdm_fn;
    saveAs_fn = params.saveAs_fn;

    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    
    % ------------------------------------------
    % Set up functional filenames for processing
    % ------------------------------------------

    % Could be a single filename, or cell array of filenames
    if iscell(functional_fn)
        N_runs = numel(functional_fn);
        %func_files = functional_fn;
        if ~iscell(saveAs_fn) || numel(saveAs_fn)~=N_runs
            % Add error message TODO
            disp('Error: saveAs_fn should have the same data type and number of indices as functional_fn')
        end
        
    elseif isstring(functional_fn) || ischar(functional_fn)
        fn_charstring = true;
        N_runs = 1;
        func_files = {functional_fn};
        saveAs_fn = {saveAs_fn}; % TODO first test if the save fn is also a single string/chararray
    else
        % Add error message TODO
        disp('Error: functional file input is not cell array, nor char/string')
    end
    
    for r = 1:N_runs
        vbdm_apply = struct;
        vbdm_apply.matlabbatch{1}.spm.tools.fieldmap.applyvdm.data.scans(1) = functional_fn(r);
        vbdm_apply.matlabbatch{1}.spm.tools.fieldmap.applyvdm.data.vdmfile(1) = vdm_fn(r);
        vbdm_apply.matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.pedir = 2;
        vbdm_apply.matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.which = [2 1];
        vbdm_apply.matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.rinterp = 4;
        vbdm_apply.matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.wrap = [0 0 0];
        vbdm_apply.matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.mask = 1;
        vbdm_apply.matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.prefix = 'u';

        % -------------------
        % Run SPM12 job
        % -------------------
        spm_jobman('run', vbdm_apply.matlabbatch);
    end



    % -------------
    % Format output
    % -------------
    output = struct;
    for r = 1:N_runs
        [d, f, e] = fileparts(functional_fn{r});
        temp_functional_fn = fullfile(d, ['u' f e]);
        [status, msg, msgID] = movefile(temp_functional_fn, saveAs_fn{r});
        
    end    