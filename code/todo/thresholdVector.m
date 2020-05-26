function output = thresholdVector(vec, lower, upper, type)

% Type = 0 ==> set thresholded value to zero
% Type = 1 ==> set thresholded value to upper or lower threshold

if type == 0
    vec(vec <= lower) = 0;
    vec(vec >= upper) = 0;
else
    vec(vec <= lower) = lower;
    vec(vec >= upper) = upper;
end

output = vec;

end
