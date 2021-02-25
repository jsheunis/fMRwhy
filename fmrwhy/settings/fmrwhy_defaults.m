function options = fmrwhy_defaults()
% Outputs the ``options`` structure with the fMRwhy toolbox root directory location as a single field

current_path = mfilename('fullpath');
ind = strfind(current_path,'fMRwhy');
options.fmrwhy_dir = current_path(1:ind+5);