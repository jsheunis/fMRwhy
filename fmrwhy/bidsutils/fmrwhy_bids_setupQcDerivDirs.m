function options = fmrwhy_bids_setupQcDerivDirs(bids_dir, options)
% Creates (or checks) the correct BIDS derivatives directory structure for the Quality Control workflow (:func:`fmrwhy_workflow_qc`).
%
% :param bids_dir: Directory location of the BIDS dataset for which the derivatives are created.
% :type bids_dir: character array
% :param options: Empty or existing `options` structure
% :type options: struct
% :returns: ``options`` - updated structure with directory locations of located dependencies

if nargin < 2
    options = struct;
end

% Derivatives directory structure
options.deriv_dir = fullfile(bids_dir, 'derivatives');
options.preproc_dir = fullfile(options.deriv_dir, 'fmrwhy-preproc');
options.qc_dir = fullfile(options.deriv_dir, 'fmrwhy-qc');
% If the fMRwhy derivative directories do not exist, create them
dirs = {options.preproc_dir, options.qc_dir};
for i = 1:numel(dirs)
    if ~exist(dirs{i}, 'dir')
        mkdir(dirs{i});
    end
end