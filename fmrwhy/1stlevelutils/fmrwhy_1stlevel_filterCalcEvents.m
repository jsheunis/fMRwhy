function cond = fmrwhy_1stlevel_filterCalcEvents(event_fn, condition_names, cond_params)
    % 

    if ~iscell(condition_names)
        error('Condition names should be a cell array')
    end

    T = readtable(event_fn, 'FileType', 'text', 'delimiter', '\t');

    for i = 1:numel(condition_names)
        cond_name = condition_names{i};

        if ~isfield(cond_params,cond_name)
            disp(['WARNING: No filtering or calculation parameters specified for condition name: ' cond_name]);
            continue;
        end

        Nc = numel(cond_params.(cond_name));
        if Nc > 1
            filtered_T = [];
            calculated_T = [];
            for c = 1:Nc
                filter_params = cond_params.(cond_name)(c).filter;
                calc_params = cond_params.(cond_name)(c).calc;
                filt_t = fmrwhy_util_tableFilter(T, filter_params);
                calc_t = fmrwhy_util_tableCalc(filt_t, calc_params);
                filtered_T = [filtered_T; filt_t];
                calculated_T = [calculated_T; calc_t];
            end 
        else
            filter_params = cond_params.(cond_name).filter;
            calc_params = cond_params.(cond_name).calc;

            filtered_T = fmrwhy_util_tableFilter(T, filter_params);
            calculated_T = fmrwhy_util_tableCalc(filtered_T, calc_params);
        end

        if isempty(calculated_T.onset)
            disp(['Empty duration column: ' cond_name])
            continue
        else
            cond(i).name = cond_name;
            cond(i).onset = calculated_T.onset;
            cond(i).duration = calculated_T.duration;
            cond(i).tmod = 0;
            cond(i).pmod = {''};
            cond(i).orth = 1;
        end
    end


