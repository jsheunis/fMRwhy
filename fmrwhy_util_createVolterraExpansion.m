function out_data = fmrwhy_util_createVolterraExpansion(data, include_data)

% Matlab script derived from John' Gems:
% https://www2.warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/scripts/spm/mp_diffpow24.sh

% Quote: "Creates file with 24 columns; the first 6 are the motion parameters, the
% next 6 are the square of the motion parameters, the next 6 are the
% temporal difference of motion parameters, and the next 6 are the square
% of the differenced values. This is useful for accounting for 'spin
% history' effects, and variation not otherwise accounted for by motion
% correction."

% Expects Nx6 matrix
% out_data: R = [(MP), squares, derivatives, squaresOfDerivatives]

[Ni, Nj] = size(data);
if Nj ~= 6
    disp(['Number of movement parameters incorrect, or matrix transposed - ' num2str(Ni) ' detected instead of 6.']);
end

out_data = zeros(Ni, Nj*4);
out_data(:, 1:Nj) = data;
out_data(:, (Nj+1):2*Nj) = data.^2;
temp_diff = [zeros(1, Nj); diff(data)];
out_data(:, (2*Nj+1):3*Nj) = temp_diff;
out_data(:, (3*Nj+1):4*Nj) = temp_diff.^2;

if ~include_data
    [Mi, Mj] = size(out_data);
    out_data = out_data(:, (Nj+1):Mj);
end





