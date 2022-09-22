function fmrwhy_batch_VDMcalc(functional_fn, saveAs_fn, phase1_fmaps, magnitude1_fmaps, phase2_fmaps, magnitude2_fmaps, varargin)
    
    % -------------
    % Parse inputs
    % -------------
    validChar = @(x) ischar(x);
    validType = @(x) any(validatestring(x, filetypes));

%     eOptions = {'quality', 'sep', 'fwhm', 'rtm', 'einterp', 'ewrap', 'weight'};
%     eDefaults = {0.9, 4, 5, 0, 2, [0 0 0], ''};
%     rOptions = {'which', 'rinterp', 'rwrap', 'mask', 'prefix'};
%     rDefaults = {[1 0], 4, [0 0 0], 1, 'r'};
% 
    p = inputParser;
    addRequired(p, 'functional_fn');
    addRequired(p, 'saveAs_fn');
    addRequired(p, 'phase1_fmaps');
    addRequired(p, 'magnitude1_fmaps');
    addRequired(p, 'phase2_fmaps');
    addRequired(p, 'magnitude2_fmaps');
%     for i = 1:numel(eOptions)
%         addParameter(p, eOptions{i}, eDefaults{i});
%     end
%     for i = 1:numel(rOptions)
%         addParameter(p, rOptions{i}, rDefaults{i});
%     end
    parse(p, functional_fn, saveAs_fn, phase1_fmaps, magnitude1_fmaps, phase2_fmaps, magnitude2_fmaps, varargin{:});
    params = p.Results;
    functional_fn = params.functional_fn;
    phase1_fmaps = params.phase1_fmaps;
    phase2_fmaps = params.phase2_fmaps;
    magnitude1_fmaps = params.magnitude1_fmaps;
    magnitude2_fmaps = params.magnitude2_fmaps;
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
        vbdm_calc = struct;
        
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortphase = phase1_fmaps(r);
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortmag = magnitude1_fmaps(r);
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longphase = phase2_fmaps(r);
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longmag = magnitude2_fmaps(r);
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [4.6 6.9];
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = 14.256;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {'/Users/lhellr/Documents/MATLAB/AdditionalStuff/spm12/toolbox/FieldMap/T1.nii'};
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = functional_fn(r);
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 1;
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
        vbdm_calc.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 1;
        
        % -------------------
        % Run SPM12 job
        % -------------------
        spm_jobman('run', vbdm_calc.matlabbatch);
    end

%###############



    
    % Data
    %vdm5_scsub-30823_ses-2_run-1_task-NFTask_phase1
    % prefix = 'vdm5_sc'
    % appended is phase1_fmaps{r}

    % -------------
    % Format output
    % -------------
    output = struct;
    for r = 1:N_runs
        [d, f, e] = fileparts(phase1_fmaps{r});
        temp_vdm_fn = fullfile(d, ['vdm5_sc' f e]);
        [status, msg, msgID] = movefile(temp_vdm_fn, saveAs_fn{r});
        
    end