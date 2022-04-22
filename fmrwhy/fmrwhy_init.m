function fmrwhy_init(dependency_dir)
    % add fmrwhy and subdirectories to MATLAB PATH
    current_path = mfilename('fullpath');
    ind = strfind(current_path, 'fMRwhy');
    fmrwhy_dir = [current_path(1:ind + 5) filesep 'fmrwhy'];
    addpath(genpath(fmrwhy_dir));

    if nargin<1
        dependency_dir = [current_path(1:ind + 5) filesep 'dependencies'];

    addpath([dependency_dir filesep 'spm12']);
    addpath(genpath([dependency_dir filesep 'bids-matlab']));
    addpath(genpath([dependency_dir filesep 'dicm2nii']));
    addpath(genpath([dependency_dir filesep 'tapas' filesep 'physio' filesep 'code']));
    tapas_physio_init()
    addpath(genpath([dependency_dir filesep 'RainCloudPlots']));

    % check if dependencies are installed and throw warnings if not
    fmrwhy_util_checkDependencies([]);
