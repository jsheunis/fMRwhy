function fmrwhy_batch_VBDMcalcApply(rfunctional_fn, template_fn, saveAs_fn, functional_fns, phase1_fmaps, magnitude1_fmaps, phase2_fmaps, magnitude2_fmaps, varargin)
    
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
    addRequired(p, 'rfunctional_fn');
    addRequired(p, 'template_fn');
    addRequired(p, 'saveAs_fn');
    addRequired(p, 'functional_fns');
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
    parse(p, rfunctional_fn, template_fn, saveAs_fn, functional_fns, phase1_fmaps, magnitude1_fmaps, phase2_fmaps, magnitude2_fmaps, varargin{:});
    params = p.Results;
    rfunctional_fn = params.functional_fn;
    template_fn = params.template_fn;
    functional_fns = params.functional_fns;
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
    if iscell(rfunctional_fn)
        N_runs = numel(rfunctional_fn);
        func_files = rfunctional_fn;
        if ~iscell(saveAs_fn) || numel(saveAs_fn)~=N_runs
            % Add error message TODO
            disp('Error: saveAs_fn should have the same data type and number of indices as functional_fn')
        end
        
    elseif isstring(rfunctional_fn) || ischar(rfunctional_fn)
        fn_charstring = true;
        N_runs = 1;
        func_files = {rfunctional_fn};
        saveAs_fn = {saveAs_fn}; % TODO first test if the save fn is also a single string/chararray
    else
        % Add error message TODO
        disp('Error: functional file input is not cell array, nor char/string')
    end
    
    data = cell(1,N_runs);
    temp_functional_fn = cell(1,N_runs);
    for r = 1:N_runs
        vbdm_calc_apply = struct;
        
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortphase = phase1_fmaps{r};
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortmag = magnitude1_fmaps{r};
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longphase = phase2{r};
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longmag = magnitude2{r};
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [4.6 6.9];
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = 14.256;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {'/Users/lhellr/Documents/MATLAB/AdditionalStuff/spm12/toolbox/FieldMap/T1.nii'};
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = functional_fns{r};
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 1;
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
        vbdm_calc_apply.matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 1;
        
        
        vbdm_calc_apply.matlabbatch{2}.spm.tools.fieldmap.applyvdm.data.scans(1) = rfunctional_fn{r};
        %vdm5_scsub-30823_ses-2_run-1_task-NFTask_phase1.nii
        vbdm_calc_apply.matlabbatch{2}.spm.tools.fieldmap.applyvdm.data.vdmfile(1) = cfg_dep('Calculate VDM: Voxel displacement map (Subj 1, Session 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','vdmfile', '{}',{1}));
        vbdm_calc_apply.matlabbatch{2}.spm.tools.fieldmap.applyvdm.roptions.pedir = 2;
        vbdm_calc_apply.matlabbatch{2}.spm.tools.fieldmap.applyvdm.roptions.which = [2 1];
        vbdm_calc_apply.matlabbatch{2}.spm.tools.fieldmap.applyvdm.roptions.rinterp = 4;
        vbdm_calc_apply.matlabbatch{2}.spm.tools.fieldmap.applyvdm.roptions.wrap = [0 0 0];
        vbdm_calc_apply.matlabbatch{2}.spm.tools.fieldmap.applyvdm.roptions.mask = 1;
        vbdm_calc_apply.matlabbatch{2}.spm.tools.fieldmap.applyvdm.roptions.prefix = 'u';
        
        % -------------------
        % Run SPM12 job
        % -------------------
        spm_jobman('run', vbdm_calc_apply.matlabbatch);
    end

%###############



    
    % Data
    

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