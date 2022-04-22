function calculated_T = fmrwhy_util_tableCalc(T, calc_params)
    % Description
    %
    % :param T: full table with content and column headers
    % :type T: Table
    % :param calc_params: full table with content and column headers
    % :type calc_params: Cell array
    % :returns: calculated_T - the output table with calculated columns

    calculated_T = T;
    for i = 1:numel(calc_params)
        [col_name, col_vals] = columnCalc(calculated_T, calc_params{i});
        calculated_T.(col_name) = col_vals;
    end

function [col_name, col_vals] = columnCalc(Tb, params)

    p1 = params{1};
    p2 = params{2};

    if iscell(p2)
        [~, p2] = columnCalc(Tb, p2);
    end

    if ischar(p2)
        p2 = Tb.(p2);
    end
    
    if ischar(p1)
        col_name = p1(1:end-1);
        math_op = p1(end);
        p1num = getNum(Tb.(col_name));
        p2num = getNum(p2);
        if numel(p2num) == 1
            p2mat = p2num*ones(numel(p1num), 1);
        else
            p2mat = p2num;
        end

        switch math_op
            case '+'
                col_vals = p1num + p2mat;
            case '-'
                col_vals = p1num - p2mat;
            case '*'
                col_vals = p1num .* p2mat;
            case '/'
                col_vals = p1num ./ p2mat;
            case '='
                col_vals = p2mat;
            otherwise
                error('Last character of first parameter of a table column calculation should be a mathematical operation')
        end
    else
        % p1{}
        error('First parameter of a table column calculation should be a char')
    end


function nm = getNum(strnm)
    if isnumeric(strnm)
        nm = strnm;
    else
        nm =  str2num(char(strnm));
    end