function R = createMotionRegressors(MP)

% Matlab script derived from John' Gems: https://www2.warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/scripts/spm/mp_diffpow24.sh

% "Creates file with 24 columns; the first 6 are the motion parameters, the
% next 6 are the square of the motion parameters, the next 6 are the
% temporal difference of motion parameters, and the next 6 are the square
% of the differenced values. This is useful for accounting for 'spin
% history' effects, and variation not otherwise accounted for by motion
% correction."
[Ni, Nj] = size(MP);
if Nj ~= 6
    disp(['Number of movement parameters incorrect - ' num2str(Ni) ' detected.']);
end

R = zeros(Ni, Nj*4);
R(:, 1:Nj) = MP;
R(:, (Nj+1):2*Nj) = MP.^2;
temp_diff = [zeros(1, Nj); diff(MP)];
R(:, (2*Nj+1):3*Nj) = temp_diff;
R(:, (3*Nj+1):4*Nj) = temp_diff.^2;







