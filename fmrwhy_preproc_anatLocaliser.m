% FUNCTION:
%--------------------------------------------------------------------------

% Copyright statement....

%--------------------------------------------------------------------------
% DEFINITION
%--------------------------------------------------------------------------

%...

% INPUT:
% funcional0_fn     - filename of initial pre-real-time 3D functional volume template
% structural_fn     - filename of T1-weighted structural volume
% spm_dir           - SPM12 directory
%
% OUTPUT:
% output            - structure with filenames and data

%--------------------------------------------------------------------------
% STEPS

% Assumes basic structura-functional preproc has been done (coreg, etc), i.e.
% that the subject to MNI space mapping exists in the preproc_data stucture.
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------

function output = rtme_preproc_anatLocaliser(sub, defaults)

% Load required defaults
spm_dir = defaults.spm_dir;
preproc_dir = defaults.preproc_dir;
template_run = defaults.template_run;
template_task = defaults.template_task;
template_echo = defaults.template_echo;


% Grab files for preprocessing
functional0_fn = fullfile(preproc_dir, sub, 'func', [sub '_task-' template_task '_run-' template_run '_echo' template_echo 'bold.nii,1']);
inverse_transformation = fullfile(preproc_dir, sub, 'anat', ['iy_' sub '_T1w.nii']);
ROI_fns = defaults.ROI_fns;

% Structure to save outputs
output = struct;

% STEP 1 -- Warp MNI space rois to functional space
disp('1 - Warp MNI space ROIs to functional space...');
spm('defaults','fmri');
spm_jobman('initcfg');
normalize_write = struct;
% Deformation field
normalize_write.matlabbatch{1}.spm.spatial.normalise.write.subj.def = {inverse_transformation};
% Data
normalize_write.matlabbatch{1}.spm.spatial.normalise.write.subj.resample = ROI_fns';
% Write options
normalize_write.matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78,-112,-70;78,76,85];
normalize_write.matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2,2,2];
normalize_write.matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
normalize_write.matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
% Run
spm_jobman('run', normalize_write.matlabbatch);
% Save outputs
for roi = 1:numel(ROI_fns)
    [droi, fnroi, extroi] = fileparts(ROI_fns{roi});
    output.wROI_fns{roi} = [droi filesep 'w' fnroi extroi];
end
disp('done')

% STEP 2 -- Reslice warped ROIs to functional space grid
disp('2 - Reslice warped ROIs to functional space grid...');
spm('defaults','fmri');
spm_jobman('initcfg');
reslice = struct;
% Ref
reslice.matlabbatch{1}.spm.spatial.coreg.write.ref = {functional0_fn};
% Source
source_fns = output.wROI_fns;
reslice.matlabbatch{1}.spm.spatial.coreg.write.source = source_fns';
% Roptions
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
% Run
spm_jobman('run',reslice.matlabbatch);
% Save outputs
for roi = 1:numel(ROI_fns)
    [droi, fnroi, extroi] = fileparts(ROI_fns{roi});
    output.rwROI_fns{roi} = [droi filesep 'rw' fnroi extroi];
end
disp('done')