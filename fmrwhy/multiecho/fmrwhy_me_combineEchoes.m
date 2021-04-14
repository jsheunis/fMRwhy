function combined_data = fmrwhy_me_combineEchoes(func_data, TE, mask, method, weight_data)
    % function combined_img2D = fmrwhy_me_combineEchoes(TE, method, weight_img, F, I_mask, N)

    % Function to combine multi-echo timeseries (or single volume) data into a single combined timeseries (or volume)
    % Various combination methods are available. In all cases, a weighted average is computed. This is equivalent to
    % a weighted summation if weights are normalised.

    % INPUTS:
    % func_data - 5D [4D x E] timeseries or 4D [3D x E] volume data to be combined_data
    % TE        - vector of echo times (ms) in order of acquisition, e.g. [14 28 42]
    % mask      -
    % method    - method used for combination of echoes:
    %           1 = T2*-weighted average, using T2* map provided in the "weight_data" parameter (Posse)
    %           2 = tSNR-weighted average, using tSNR map provided in the "weight_data" parameter (PAID, Poser)
    %           3 = TE-weighted average, using echo times in TE (ignores "weight_data" parameter)
    %           4 = User value-weighted average, using weights in TE (ignores "weight_data" parameter)
    %
    % THEORY:
    % Weighted average = {sum(data x weights)} / {sum(weights)}
    % If data = X and weights = W:
    % Weighted average = (x1w1 + x2w2 + ... + xnwn)/(w1+w2+...+wn)

    sz = size(func_data);
    Ne = numel(TE);

    % Check if dimensions agree for number of echos in data and TE vector
    if Ne ~= sz(end)
        disp('error '); % TODO: raise error
        return
    end

    % Check validity of weight_data based on combination method
    switch method
        case 1
            % Posse T2*-weighted
            % weight_data should be a single 3D image same size as first three dimensions of func_data
            % isequal(size(A), size(B)) || (isvector(A) && isvector(B) && numel(A) == numel(B))
        case 2
            % Poser PAID tSNR-weighted
            % weight_data should be a 3D x E matrix of tSNR images, where E equals the amount of echoes
        case 3
            % TE-weighted

        case 4
            % User specified vector of weights, same size as TE

        otherwise
            disp('error '); % TODO: raise error
    end

    % Method 1: Per-voxel T2star-weighted combination. This 3D T2* image could
    % be calculated in various ways (e.g. T2star and S0 estimated per voxel
    % from time-series average of multiple echoes; or T2star and S0 estimated
    % in real-time). See Posse et al and Poser et al for details.
    % S = sum[S(TEn).w(TEn)n]
    % w(TEn)n = [TEn.exp(-TEn/T2*)]/[TEn.exp(?TEn/T2*)]
    % thus, S = sum[S(TEn)].[TEn.exp(?TEn/T2*)]/[TEn.exp(-TEn/T2*)]

    % Method 2: Per-voxel tSNR-weighted combination. This requires a 3D tSNR
    % image per echo. See poser et al for details.
    % S = sum[S(TEn).w(TEn)n]
    % w(TEn)n = [tSNR.TEn]/sum[tSNR.TEn]

    % Method 3: TE-weighted combination. This requires a 3D tSNR
    % image per echo. See poser et al for details.
    % S = sum[S(TEn).w(TEn)n]
    % w(TEn)n = [tSNR.TEn]/sum[tSNR.TEn]

    data = {};
    weights = {};
    weighted_data = {};
    % sum_weights = zeros(sz(1:3));
    % sum_weighted_data = zeros(sz(1:4));
    sum_weights = 0;
    sum_weighted_data = 0;

    for e = 1:Ne
        if numel(sz) == 4
            data{e} = func_data(:, :, :, e);
        else % numel(sz) = 5
            data{e} = func_data(:, :, :, :, e);
        end

        switch method
            case 1
                % Posse T2*-weighted

                weights{e} = TE(e) .* exp(-TE(e) ./ weight_data);
            case 2
                % Poser PAID tSNR-weighted
                weights{e} = TE(e) .* weight_data(:, :, :, e);
            case 3
                % TE-weighted
                weights{e} = TE(e);
            case 4
                %
                weights{e} = weight_data(e);
            otherwise
                disp('Unknown combinaion method'); % TODO: raise error
        end
        % x1w1, x2w2, ..., xnwn
        weighted_data{e} = data{e} .* weights{e};
        % w1 + w2 + wn
        sum_weights = sum_weights + weights{e};
        % x1w1 + x2w2 + ... + xnwn
        sum_weighted_data = sum_weighted_data + weighted_data{e};
    end
    % (x1w1 + x2w2 + ... + xnwn)/(w1+w2+...+wn)
    combined_data = sum_weighted_data ./ sum_weights; % a 4D timeseries [3D x t] or a 3D image

    % Mask data

    %
    % Ne = numel(TE);
    % if numel(F) ~= Ne
    %    disp('The the number of echoes in the functional image (F) does not match the number of echo times')
    %    return;
    % end
    %
    %
    %% Nargs = nargin;
    %
    %% if Nargs-4 ~= Ne
    %%     disp('Number of image inputs do not match echo time vector size')
    %% end
    %
    % combined_img2D = zeros(N(1)*N(2)*N(3), 1);
    %% sum_weights = combined_img;
    %% weights_denom = cell(Ne,1);
    %% weights = cell(Ne,1);
    %
    % switch method
    %    case 1
    %        T2star_img = weight_img;
    %        T2star_2D = reshape(T2star_img, N(1)*N(2)*N(3),1);
    %        weights_denom = TE.*exp(-TE./T2star_2D(I_mask)); % matrix = Nmask * Ne (each column represents denominator for specific echo)
    %        weights_sum = sum(weights_denom, 2); % matrix = Nmask * 1 (column represents sum of denominators for all echoes)
    %        weights = weights_denom./weights_sum; % matrix = Nmask * Ne (each column represents weighting vector for specific echo volume)
    %
    %        f = [];
    %        for e = 1:Ne
    %            fe = reshape(F{e}, N(1)*N(2)*N(3),1);
    %            f = cat(2, f, fe(I_mask));
    %        end
    %
    %        F_weighted = f.*weights;
    %        combined_masked = sum(F_weighted, 2);
    %        combined_img2D(I_mask) = combined_masked;
    %
    %    case 2
    %        if numel(weight_img) ~= Ne
    %            disp('For method 2 (tSNR-weighted combination), the weight image should have the same dimension as the number of echoes')
    %            return;
    %        end
    %
    %        f = [];
    %        tSNR = [];
    %        for e = 1:Ne
    %            fe = reshape(F{e}, N(1)*N(2)*N(3),1);
    %            f = cat(2, f, fe(I_mask));
    %            tSNRe = reshape(weight_img{e}, N(1)*N(2)*N(3),1);
    %            tSNR = cat(2, tSNR, tSNRe(I_mask));
    %        end
    %
    %        % w(TEn)n = [tSNR.TEn]/sum[tSNR.TEn]
    %        weights_denom = TE.*tSNR; % matrix = Nmask * Ne (each column represents denominator for specific echo)
    %        weights_sum = sum(weights_denom, 2); % matrix = Nmask * 1 (column represents sum of denominators for all echoes)
    %        weights = weights_denom./weights_sum; % matrix = Nmask * Ne (each column represents weighting vector for specific echo volume)
    %
    %        F_weighted = f.*weights;
    %        combined_masked = sum(F_weighted, 2);
    %        combined_img2D(I_mask) = combined_masked;
    %
    %    otherwise
    %        disp('No implementation for this method')
    % end
