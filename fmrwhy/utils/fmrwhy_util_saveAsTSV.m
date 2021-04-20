function tsv_fn = fmrwhy_util_saveAsTSV(txt_fn, col_names)

    if ~isempty(txt_fn)
        [d, fn, ext] = fileparts(txt_fn);
        tsv_fn = fullfile(d, [fn '.tsv']);
        temp_txt_fn = fullfile(d, [fn '_temp.txt']);

        data = load(txt_fn);
        data_table = array2table(data, 'VariableNames', col_names);
        writetable(data_table, temp_txt_fn, 'Delimiter', '\t');
        [status, msg, msgID] = movefile(temp_txt_fn, tsv_fn);
    else
    end
