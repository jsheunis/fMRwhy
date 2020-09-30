function out_data = fmrwhy_util_demean(data)
% Removes timeseries mean from data
%
% INPUT:
% data    - R x C matrix; R = rows = time points;
%           C = columns = variables/parameters
%
% OUTPUT:
% out_data

% Define variables
[r, c] = size(data);
% Remove mean from data
data = data - repmat(nanmean(data, 1), r, 1);
out_data = data;