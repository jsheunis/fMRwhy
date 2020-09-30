function fmrwhy_util_checkDependencies()

spm_url = 'https://github.com/spm/spm12/releases/tag/r7771';
bidsmatlab_url = 'https://github.com/bids-standard/bids-matlab';
dicm2nii_url = 'https://github.com/jsheunis/dicm2nii/releases/tag/v0.2';
tapasphysio_url = 'https://github.com/translationalneuromodeling/tapas/releases/tag/v4.0.0';
raincoud_url = 'https://github.com/RainCloudPlots/RainCloudPlots/releases/tag/v1.1';

% Check and add SPM path
spm_installed = exist('spm');
if ~spm_installed
    msg = 'SPM12 cannot be found on your MATLABPATH. Please add the top-level SPM12 directory to the path using e.g. `addpath(/path/to/your/spm12/directory)`';
    errorStruct.identifier = 'checkDependencies:missingDependency';
    errorStruct.message = sprintf('%s \n%s  \n%s', ...
        'SPM12 cannot be found on your MATLABPATH.', ...
        'Please add the top-level SPM12 directory to the path using e.g. `addpath(/path/to/your/spm12/directory)`.', ...
        'You can download the required release of SP12 here: %s', spm_url);
    error(errorStruct);
else
    pathSpm = spm('Dir');
    addpath(pathSpm)
    spm_check_installation();
end

% Check bids-matlab
% No unique function names to check for using which, so using `exist('bids-matlab')` for now
bm_installed = exist('bids-matlab');
if ~bm_installed
    createErrorMsg('bids-matlab', bidsmatlab_url);
else
    bm_what = what('bids-matlab');
    addpath(genpath(bm_what.path))
end

% Check dicm2nii
dicm2nii_installed = which('dicm2nii');
if isempty(dicm2nii_installed)
    createErrorMsg('dicm2nii', dicm2nii_url)
else
    [dicm2nii_dir, f, e] = fileparts(dicm2nii_installed);
    addpath(genpath(dicm2nii_dir))
end

% Check TAPAS physio
tapasphysio_installed = which('tapas_init');
if isempty(tapasphysio_installed)
    createErrorMsg('TAPAS', tapasphysio_url)
else
    [tapasphysio_dir, f, e] = fileparts(tapasphysio_installed);
    addpath(genpath(tapasphysio_dir))
end

% Check raincloudplots
raincloud_installed = which('raincloud_plot');
if isempty(raincloud_installed)
    createErrorMsg('Raincloud plots', raincloud_url)
else
    [raincloud_dir, f, e] = fileparts(raincloud_installed);
    addpath(genpath(raincloud_dir))
end

%-----------------
% Helper functions
%-----------------

 function createErrorMsg(name, url)
    errorStruct.identifier = 'checkDependencies:missingDependency';
    errorStruct.message = sprintf('\n%s \n\n%s  \n\n%s', ...
        [name ' cannot be found on your MATLABPATH.'], ...
        ['Please add the top-level ' name ' directory and all subdirectories to the path using e.g. `addpath(genpath(/path/to/your/' name '/directory))`.'], ...
        ['You can download the required release of ' name ' here: ' url]);
    error(errorStruct);
 