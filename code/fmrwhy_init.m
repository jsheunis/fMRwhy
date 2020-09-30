% add fmrwhy and subdirectories to MATLAB PATH
current_path = mfilename('fullpath');
ind = strfind(current_path,'fMRwhy');
fmrwhy_dir = current_path(1:ind+5);
addpath(genpath(fmrwhy_dir))

% check if dependencies are installed and throw warnings if not
fmrwhy_util_checkDependencies();