function [filename, filepath] = fmrwhy_bids_constructFilename(filetype, varargin)
    % BIDS 1.4.0 - This does not yet take into account whether an entity is OPTIONAL or REQUIRED for a specific file modality/type (e.g. 'func', 'anat', 'meg', etc).
    % It assumes all parameters are optional and just constructs the filename given the input also, file extension is passed as parameter ==> should be automatically determined from file modality/type
    %
    % :param filetype: something x
    % :param varargin: something y
    % :returns: filename

    filetypes = {'anat', 'func', 'fmap'}; % currently supported in fMRwhy
    descriptions = {'Subject', 'Session', 'Task', 'Acquisition', 'Contrast Enhancing Agent', 'Reconstruction', 'Phase-Encoding Direction', 'Run', 'Corresponding modality', 'Echo', 'Recording', 'Processed (on device)', 'Space', 'Split', 'Description', 'File extension'};
    entities = {'sub', 'ses', 'task', 'acq', 'ce', 'rec', 'dir', 'run', 'mod', 'echo', 'recording', 'proc', 'space', 'split', 'desc', 'ext'};
    formats = {'label', 'label', 'label', 'label', 'label', 'label', 'label', 'index', 'label', 'index', 'label', 'label', 'label', 'index', 'label', 'label'};
    % {'sub', 'ses', 'acq', 'ce', 'rec', 'fa', 'echo', 'inv', 'run'} % SPM12 list seems to be pre-v1.4.0

    validChar = @(x) ischar(x);
    validType = @(x) any(validatestring(x, filetypes));

    p = inputParser;
    addRequired(p, 'filetype', validType);
    for i = 1:numel(entities)
        addParameter(p, entities{i}, '', validChar);
    end
    parse(p, filetype, varargin{:});

    params = p.Results;
    filename = '';
    for i = 1:numel(entities)
        if ~isempty(params.(entities{i}))
            if i == 1
                newStr = [entities{i} '-' params.(entities{i})];
            elseif i == numel(entities) % last entity = extension
                newStr = [params.(entities{i})];
            else
                newStr = ['_' entities{i} '-' params.(entities{i})];
            end
            filename = [filename newStr];
        end
    end

    filepath = ['sub-' params.sub];
    if ~isempty(params.ses)
        filepath = [filepath filesep 'ses-' params.ses];
    end
    if ~isempty(params.filetype)
        filepath = [filepath filesep params.filetype];
    end
