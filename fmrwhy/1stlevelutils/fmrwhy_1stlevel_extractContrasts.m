function consess = fmrwhy_1stlevel_extractContrasts(SPM, contrast_params, runs)
    % Vertically concatenate the columns present in a set of events files that belong to 
    % a set of functional files, typically consecutive runs of the same task.
    %
    % :param SPM: 
    % :type SPM: struct
    % :param contrast_params: 
    % :type contrast_params: struct
    % :param runs: 
    % :type runs: cell array
    % :returns: consess - cell array

    [Ntt, Nregr] = size(SPM.xX.X);
    all_regressor_names = SPM.xX.name;
    bf_extension = '*bf(1)';
    str_1 = 'Sn(';
    str_2 = ') ';
    consess = {};
    for c=1:numel(contrast_params)
        % Match name and add + subtract conditions
        matched_name_add = zeros(1, Nregr);
        matched_name_subtract = zeros(1, Nregr);

        for r = 1:numel(runs)
            rn = runs{r};
            % add conditions
            if isfield(contrast_params(c),'add') && ~isempty(contrast_params(c).add)
                for a = 1:numel(contrast_params(c).add)
                    str_to_match_add = [str_1 rn str_2 contrast_params(c).add{a} bf_extension];
                    matched_name_add = matched_name_add | strcmp(all_regressor_names, str_to_match_add);
                end
            end
            % subtract conditions
            if isfield(contrast_params(c),'subtract') && ~isempty(contrast_params(c).subtract)
                for s = 1:numel(contrast_params(c).subtract)
                    str_to_match_subtract = [str_1 rn str_2 contrast_params(c).subtract{s} bf_extension];
                    matched_name_subtract = matched_name_subtract | strcmp(all_regressor_names, str_to_match_subtract);
                end
            end
        end
        % Create enw contrast
        consess{c}.tcon.name = contrast_params(c).name;
        consess{c}.tcon.weights = matched_name_add - matched_name_subtract;
        consess{c}.tcon.sessrep = 'none';
    end