% fmrwhy_script_calculatePSCfromBeta

beta_effect = '/Volumes/TSM/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-motor_run-1_echo-2/beta_0001.nii';
beta_constant = '/Volumes/TSM/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-motor_run-1_echo-2/beta_0019.nii';
spmmat = '/Volumes/TSM/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-motor_run-1_echo-2/SPM.mat';

spm = load(spmmat);
ntime = 20*(1/spm.SPM.xBF.dt);
reference_block =  conv(ones(1,ntime),SPM.xBF.bf(:,1))'
SF = max(reference_block);


nii = nii_tool('load', beta_effect);
beta_vals = double(nii.img);
nii = nii_tool('load', beta_constant);
const_vals = double(nii.img);

PSC_vals = beta_vals * SF ./ const_vals * 100;

new_fn = '/Volumes/TSM/NEUFEPME_data_BIDS/derivatives/fmrwhy-stats/sub-001/task-motor_run-1_echo-2/PSC_0001.nii';
no_scaling = 1;
fmrwhy_util_saveNifti(new_fn, PSC_vals, beta_effect, no_scaling)
