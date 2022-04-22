function filtered_T = fmrwhy_util_tableFilter(T, filter_params)
    % Description
    %
    % :param T: full table with content and column headers
    % :type T: Table
    % :param filter_params: full table with content and column headers
    % :type filter_params: Cell array
    % :returns: filtered_T - the filtered table

    filtered_T = T;
    for i = 1:numel(filter_params)
        filtered_T = filtered_T(strcmp(filtered_T.(filter_params{i}{1}), filter_params{i}{2})==1, :);
    end