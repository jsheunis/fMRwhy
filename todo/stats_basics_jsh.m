%% Data
x1 = randn(1,100);
x2 = randn(1,100);
%% Script
figure; plot(1:100, x1, 'bo'); hold on; plot(1:100, x2, 'ro'); hold off; legend; 
corrP = corrP_jsh(x1, x2)
[a, b] = corr(x1, x2)

%% Mean
function mu = mean_jsh(data_in)
sum_data = sum(data_in);
N = length(data_in);
mu = sum_data/N;
end

%% Variance (of sample)
function v = var_jsh(data_in)
mu = mean_jsh(data_in);
sq_diff = (data_in - mu).^2;
N = length(data_in);
v = sum(sq_diff)/(N-1);
end

%% Standard deviation
function stddev = stddev_jsh(data_in)
v = var_jsh(data_in);
stddev = sqrt(v);
end


%% Correlation
function corrP = corrP_jsh(data_in1, data_in2)

mu1 = mean_jsh(data_in1);
mu2 = mean_jsh(data_in2);
diff1 = (data_in1 - mu1);
diff2 = (data_in2 - mu2);
N1 = length(data_in1);
N2 = length(data_in2);

corrP =  sum(diff1.*diff2)/sqrt(sum(diff1.^2)*sum(diff2.^2));
end




