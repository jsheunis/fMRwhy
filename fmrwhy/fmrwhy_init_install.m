% add fmrwhy and subdirectories to MATLAB PATH
current_path = mfilename('fullpath');
ind = strfind(current_path, 'fMRwhy');
fmrwhy_dir = [current_path(1:ind + 5) filesep 'fmrwhy'];

% Dependencies
spm_url = 'https://github.com/spm/spm12/releases/tag/r7771';
bidsmatlab_url = 'https://github.com/jsheunis/bids-matlab/releases/tag/v0.0.2';
dicm2nii_url = 'https://github.com/jsheunis/dicm2nii/releases/tag/v0.2';
tapasphysio_url = 'https://github.com/translationalneuromodeling/tapas/releases/tag/v4.0.0';
raincloud_url = 'https://github.com/RainCloudPlots/RainCloudPlots/releases/tag/v1.1';


addpath(genpath(fmrwhy_dir));
dependency_dir = [current_path(1:ind + 5) filesep 'dependencies'];
addpath([dependency_dir filesep 'spm12']);
addpath(genpath([dependency_dir filesep 'bids-matlab']));
addpath(genpath([dependency_dir filesep 'dicm2nii']));
addpath(genpath([dependency_dir filesep 'tapas' filesep 'physio' filesep 'code']));
tapas_physio_init()
addpath(genpath([dependency_dir filesep 'RainCloudPlots']));

% check if dependencies are installed and throw warnings if not
fmrwhy_util_checkDependencies([]);
