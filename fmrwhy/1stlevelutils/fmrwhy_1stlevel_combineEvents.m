function combined_events_struct = fmrwhy_1stlevel_combineEvents(func_fns, event_fns, TR, saveAs_fn)
    % Vertically concatenate the columns present in a set of events files that belong to 
    % a set of functional files, typically consecutive runs of the same task.
    %
    % :param func_fns: Paths where functional files are located
    % :type func_fns: cell array
    % :param event_fns: Paths where events files are located, which should belong to the respective functional files specified in ``func_fns``
    % :type event_fns: cell array
    % :param TR: Repetition time in seconds
    % :type TR: double
    % :param saveAs_fn: Path to which the combined events information will be saved, if specified.
    % :type saveAs_fn: string or character array
    % :returns: combined_events_struct -

    if numel(func_fns) ~= numel(func_fns)
        error('The same number of functional files and event files have to be provided as arguments to fmrwhy_1stlevel_combineEvents()')
    end

    saveAs_fn_txt = '';
    if nargin == 4 && ~isempty(saveAs_fn)
        if strfind(saveAs_fn, 'txt')
            saveAs_fn_txt = saveAs_fn;
            saveAs_fn_tsv = strrep(saveAs_fn_txt, 'txt', 'tsv');
        elseif strfind(saveAs_fn, 'tsv')
            saveAs_fn_tsv = saveAs_fn;
            saveAs_fn_txt = strrep(saveAs_fn_tsv, 'tsv', 'txt');
        else
            error('saveAs_fn file extension not recognised: only txt or tsv files are allowed');
        end
    end

    end_time = 0;
    combined_events_struct = struct;

    for i = 1:numel(func_fns)
        % Get number of volumes from functional
        func_fn = func_fns{i};
        func4D_spm = spm_vol(func_fn);
        Nt = numel(func4D_spm);
        % load events into structure
        event_fn = event_fns{i};
        events_struct = tdfread(event_fn);
        % Get keys from events structure
        keys = fieldnames(events_struct);
        % For first events file, set combined_events_struct equal to the first events_struct
        % For all others, concatenate all values vertically
        % When concatenating onset times, remember to add total duration of prior scans to all elements of current event file onset times
        if i==1
            combined_events_struct = events_struct;
        else
            for j = 1:numel(keys)
                key = keys{j};
                val = events_struct.(key);
                if strfind(key, 'onset')
                    val = val + end_time
                end
                combined_events_struct.(key) = [combined_events_struct.(key); val];
            end
        end
        end_time = end_time + Nt*TR;
    end

    % save combined_events to file, if specified
    if ~isempty(saveAs_fn_txt)
        combined_events = struct2table(combined_events_struct)
        writetable(combined_events, saveAs_fn_txt, 'delimiter', '\t')
        [status,msg,msgID] = copyfile(saveAs_fn_txt, saveAs_fn_tsv)
    end
