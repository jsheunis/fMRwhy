function options = fmrwhy_defaults()

% fMRwhy toolbox root directory
current_path = mfilename('fullpath');
ind = strfind(current_path,'fMRwhy');
options.fmrwhy_dir = current_path(1:ind+5);