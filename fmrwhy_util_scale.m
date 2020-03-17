function out_data = fmrwhy_util_scale(data, lower, upper)
% Scales a data vector
%
% INPUT:
%
% OUTPUT:
% out_data
%__________________________________________________________________________
% Copyright (C) Stephan Heunis 2018

minVal = min(data);
maxVal = max(data);

out_data = (data - minVal) * (upper - lower) / (maxVal - minVal) + lower