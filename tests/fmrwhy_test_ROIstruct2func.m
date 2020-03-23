% A custom workflow that does ...


%--------------------------------------------------------------------------

% Load/create required parameters
bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';

% Setup fmrwhy bids directories on workflow level
fmrwhy_defaults_setupDerivDirs(bids_dir);

% Grab default workflow params
wf_params = fmrwhy_defaults_workflow(bids_dir);

% Loop through subjects, sessions, tasks, runs, etc
sub = '001';

% Setup fmrwhy bids directories on subject level (this copies data from bids_dir)
fmrwhy_defaults_setupSubDirs(bids_dir, sub);

% Update workflow params with subject anatomical derivative filenames
wf_params = fmrwhy_defaults_subAnat(bids_dir, sub, wf_params);

% Loop through sessions, tasks, runs, etc
ses = '';
task = 'motor';
run = '1';
echo = '2';

% Update workflow params with subject functional derivative filenames
wf_params = fmrwhy_defaults_subFunc(bids_dir, sub, ses, task, run, echo, wf_params);

%defaults.ROI_dir = '/Volumes/Stephan_WD/NEUFEPME_data_templates';
%defaults.ROI_fns = {'',
%                   {{'Left_Motor_4a_4p.nii', 'Right_Motor_4a_4p.nii'}},
%                   {'Left_Amygdala_allregions.nii', 'Right_Amygdala_allregions.nii'}} % cell array of cell arrays, in same order as tasks
%defaults.ROI_names = {'',
%                   {{'Left Motor', 'Right Motor'}},
%                   {'Left Amygdala', 'Right Amygdala'}} % cell array of cell arrays, in same order as tasks
%%

%if wf_params.has_ROIs
%end
toTransform_fns = {};
saveAs_fns = {};
count = 0;
for i = 1:numel(wf_params.roi_orig_fns)
    if ~isempty(wf_params.roi_orig_fns{i})
        for j = 1:numel(wf_params.roi_orig_fns{i})
            count = count + 1;
            toTransform_fns{count} = fullfile(wf_params.roi_orig_dir, wf_params.roi_orig_fns{i}{j});
            saveAs_fns{count} = fullfile(wf_params.anat_dir_preproc, ['sub-' sub '_space-individual_desc-' wf_params.roi_desc{i}{j} '_roi.nii']);
        end
    end
end

transformation_fn = wf_params.mni_to_indiv_fn;
template_fn = wf_params.template_fn;
%%
fmrwhy_batch_normaliseWrite(toTransform_fns, transformation_fn, template_fn, saveAs_fns)

% Reslice to functional-resolution image grid
%%
reslice_fns = {};
count = 0;
for i = 1:numel(wf_params.roi_orig_fns)
    if isempty(wf_params.roi_orig_fns{i})
        roi_fns{i} = '';
        rroi_fns{i} = '';
    else
        for j = 1:numel(wf_params.roi_orig_fns{i})
            count = count + 1;
            roi_fns{i}{j} = fullfile(wf_params.anat_dir_preproc, ['sub-' sub '_space-individual_desc-' wf_params.roi_desc{i}{j} '_roi.nii']);
            reslice_fns{count} = roi_fns{i}{j};
            rroi_fns{i}{j} = fullfile(wf_params.anat_dir_preproc, ['sub-' sub '_space-individual_desc-r' wf_params.roi_desc{i}{j} '_roi.nii']);
            saveAs_fns{count} = rroi_fns{i}{j};
        end
    end
end
%%
fmrwhy_batch_coregResl(reslice_fns, template_fn, saveAs_fns)
%%

overlay_img = fmrwhy_util_createBinaryImg(rroi_fns{2}{1}, 0); % left_motor_roi_img
template_img = spm_read_vols(spm_vol(template_fn));
output = fmrwhy_util_createOverlayMontage(template_img, overlay_img, 9, 0, 'overlayimg', 'gray', 'on', 'max');


