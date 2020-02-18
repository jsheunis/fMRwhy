dof = 10; % degrees of freedom
x_bound = 5; % T-value bounds

x = -x_bound:0.1:x_bound;
y = tpdf(x,dof);
z = normpdf(x,0,1);
figure;
plot(x,y,'-',x,z,'-.')
title(['Two-tailed probability density of Students T-distribution (dof = ' num2str(dof) ')']);
xlabel('t-value')
ylabel('Probability of observing t')



%% Data

figure;
Npoints = 1000;
trueA = 5;
trueB = 6;
Gnoise = randn(Npoints,1);
dataA = trueA + Gnoise;
Gnoise = randn(Npoints,1);
dataB = trueB + Gnoise;
subplot(211)
plot(dataA)
hold on
plot(dataB,'r')
hold off
subplot(212)
numbins = 100;
h1 = histogram(dataA, 'NumBins', numbins, 'FaceColor','b','EdgeColor','k','facealpha',0.5,'edgealpha',0.5);
hold on
h2 = histogram(dataB, 'NumBins', numbins, 'FaceColor','r','EdgeColor','k','facealpha',0.5,'edgealpha',0.5);
hold off


dataA_mean = sum(dataA)/Npoints;
dataB_mean = sum(dataB)/Npoints;
dataA_var = sum((dataA - dataA_mean).^2)/(Npoints-1); % sample variance = s^2 = sum((x - xmean)^2)/(N-1)
dataA_std = sqrt(dataA_var); % sample standard deviation = s = sqrt(var)
dataB_var = sum((dataB - dataB_mean).^2)/(Npoints-1); % sample variance = s^2 = sum((x - xmean)^2)/(N-1)
dataB_std = sqrt(dataB_var); % sample standard deviation = s = sqrt(var)


% t-test and significance testing

% We assume that the population follows a normal distribution

% Is it likely/probable that the collected sample comes from a population where
% the mean is the same as the mean of the formulated null hypothesis?

% We compare the mean of the collected sample to the null hypothesis mean

% We have to determine: is it likely that our sample mean is actually
% different from the expected null mean, or is it just a result of random
% sampling or measurement error. I.e. is the difference in the means large
% enough and is the standard deviation / error small enough so that the one
% won't be confused with the other.

% We calculate T, which is the number of standard errors that the sample
% mean is removed from the null hypothesis mean. 

% We don't have the value of the standard error (deviation), because we
% don't have the full population data. So we estimate the standard error
% from the sample(s?). Because we estimate the error, we introduce an extra
% error. Because of this we now have to assume/use a Student's T
% distribution and not the normal Z distribution anymore (more explanation needed here!)

% T = (sample mean - null mean)/(standard error);
% standard error = (sample std_dev)/sqrt(N)
sample_std = std(dataA-dataB)
se = sample_std/sqrt(Npoints)
T = (dataA_mean - dataB_mean)/se
p = 1-tcdf(T,Npoints-1) % Probability of larger t-statistic



% [h, p, ci, stats] = ttest(dataA - dataB) % paired t-test
[h, p, ci, stats] = ttest(dataA , dataB) % two-sample t-test



%% GLM sim

% From cyril pernet paper misconceptions

% 1: y = X? + ?
% with y the time series from one voxel, X the design matrix, ? the model parameters, ? the error (or residuals)
% 2: ?hat=(X_Trans.X)^(?1).X_Trans.y
% 3: ?hat^2=(ehat_Trans.ehat)/(n?rank(X))
% with ?hat the parameter estimates, ?hat^2 the variance estimate, and ehat the estimated residuals (y?X?hat)
% note that Equation 2 only applies when XTX is invertible. When XTX is rank deficient, a pseudo-inverse is used instead.


% 4: t= (c_Trans.?hat)/sqrt(?hat^2.c_Trans.(X_Trans.X)^(?1).c)
% c defined the contrast of interest, ?? are the parameter estimates, ??2 is the variance obtained from the residuals, X is the design matrix. 


