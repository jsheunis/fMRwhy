function [regressors_mat, regressors_names] = fmrwhy_1stlevel_createRegressors(confounds_fn, regressors_fn, glm_settings, run, runs)
    % :param confounds_fn: 
    % :type confounds_fn: string or character array
    % :param regressors_fn: 
    % :type regressors_fn: string or character array
    % :param options: Existing `options` structure populated with ``.firstlevel.glm_settings.(key)`` fields
    % :type options: struct
    % :returns: ``options`` - updated structure with directory locations of located dependencies

    
    % Load multiple confound regressors
    confounds_struct = tdfread(confounds_fn);
    confounds_mat = struct2array(confounds_struct);
    regressors_mat = [];
    regressors_names = {};

    % Get keys from regressor structure
    fields = fieldnames(glm_settings);
    % Process each key to see if it should be included in design matrix
    for i = 1:numel(fields)
        key = fields{i};
        val = glm_settings.(key);
        % If the value is true or larger than zero (for retroicor order),
        % parse key and include relevant data in regressor matrix
        if val
            if strfind(key, 'trans_rot')
                % Option a: any of the realignment parameters and their expansions
                if strcmp(key, 'trans_rot')
                    trans_rot_keys = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
                elseif strcmp(key, 'trans_rot_derivative1')
                    trans_rot_keys = {'trans_x_derivative1', 'trans_y_derivative1', 'trans_z_derivative1', 'rot_x_derivative1', 'rot_y_derivative1', 'rot_z_derivative1'};
                elseif strcmp(key, 'trans_rot_power2')
                    trans_rot_keys = {'trans_x_power2', 'trans_y_power2', 'trans_z_power2', 'rot_x_power2', 'rot_y_power2', 'rot_z_power2'};
                elseif strcmp(key, 'trans_rot_derivative1_power2')
                    trans_rot_keys = {'trans_x_derivative1_power2', 'trans_y_derivative1_power2', 'trans_z_derivative1_power2', 'rot_x_derivative1_power2', 'rot_y_derivative1_power2', 'rot_z_derivative1_power2'};
                else
                    % do nothing here
                end
                % Now include the applicable values
                for j = 1:numel(trans_rot_keys)
                    regressors_mat = [regressors_mat confounds_struct.(trans_rot_keys{j})];
                    regressors_names = [regressors_names {trans_rot_keys{j}}];
                end
            elseif strcmp(key, 'retroicor_c') || strcmp(key, 'retroicor_r') || strcmp(key, 'retroicor_cxr')
                % Option b: any retroicor regressors
                for k = 1:val
                    new_key = [key num2str(k)];
                    regressors_mat = [regressors_mat confounds_struct.(new_key)];
                    regressors_names = [regressors_names {new_key}];
                end
            elseif strcmp(key, 'run_number')
                continue;
            else
                % Option c: any other regressors are included as is
                regressors_mat = [regressors_mat confounds_struct.(key)];
                regressors_names = [regressors_names {key}];
            end
        end
    end

    % If glm_settings.run_number exists and is true, add run regressors
    [n, m] = size(regressors_mat);
    if any(strcmp(glm_settings,'run_number')) && glm_settings.run_number
        for i = 1:numel(runs)
            r = runs{i}
            if r == run
                add_column = ones(n,1)
            else
                add_column = zeros(n,1)
            end
            regressors_mat = [regressors_mat add_column];
            regressors_names = [regressors_names {r}];
        end
    end

    % save to files if specified
    if regressors_fn
        dlmwrite(regressors_fn, regressors_mat, 'delimiter', '\t', 'precision', '%1.7e');
        names_fn = strrep(regressors_fn, '.txt', '_names.txt');

        T = array2table(regressors_mat);
        T.Properties.VariableNames = regressors_names;
        writetable(T, names_fn, 'delimiter', '\t');

        names_fn_tsv = strrep(names_fn, '.txt', '.tsv');
        movefile(names_fn, names_fn_tsv)
    end