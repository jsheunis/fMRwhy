function options = fmrwhy_bids_getQcDerivDirs(bids_dir, options)
    % Sets the correct BIDS derivatives directory structure for the Quality Control workflow (:func:`fmrwhy_workflow_qc`).
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