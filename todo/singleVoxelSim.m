Npoints = 10000;
Gnoise = randn(Npoints,1);
subplot(321)
plot(Gnoise);
title('Gaussian noise')
subplot(322)
hist(Gnoise,20);
title('Gaussian noise histogram');
Unoise = rand(Npoints,1);
subplot(323)
plot(Unoise)
title('Uniform white noise');
subplot(324)
hist(Unoise,20);
title('Uniform white noise histogram');
% Now make auto-regressive noise
Anoise = Gnoise;
rho = 0.9;
for n=2 : length(Anoise)
 Anoise(n) = rho * Anoise(n-1) + Anoise(n);
end
subplot(325)
plot(Anoise)
title('AR(1) noise');
subplot(326)
hist(Anoise,20);
title('AR(1) noise histogram');

%%

figure(2)
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

[h, p, ci, stats] = ttest(dataA - dataB) % paired t-test
[h, p, ci, stats] = ttest(dataA , dataB) % two-sample t-test

subplot(212)
histA = hist(dataA, 20);
histB = hist(dataB, 20);
plot(histA);
hold on;
plot(histB,'r');
hold off
legend('histogram of A', 'histogram of B')
text(0.1,0.6,['degrees of freedom= ' num2str(stats.df)],'Units','Normalized');
text(0.1,0.8,['T statistic =' num2str(stats.tstat)],'Units','Normalized');
text(0.1,0.7,['Estimated Std. Dev. =', num2str(stats.sd)],'Units','Normalized');
text(0.1,0.5,['p-value =' num2str(p)],'Units','Normalized');

%%
% hrf convolved stick function - single stimulus


sr = 100; % samples per second
TR = 1/sr;
duration = 50; % seconds
hrf = spm_hrf(TR);
xBF.name = 'hrf';
xBF.dt = TR;
hrf2 = spm_get_bf(xBF);


figure;
subplot 311, plot(hrf),title('Response to a single stimulus')

% create a set of stimuli occurring at specific times
% in the time series. First make an array of zeros,
% then put ones where you want to have activity occur
stim=zeros(duration*sr,1);
stim( 5*sr:13*sr:end ) =1;
subplot 312, plot(stim),title('A set of stimuli at randomized times')
% convolve the input and the HRF and clip out the end (empty)
% first note what the convolution operation does
resp = conv(stim, hrf);
% Note the size s of the inputs and the outputs to the function:
whos resp stim hrf
% clip off the end of the response function and plot it
resp = resp(1:length(stim));
resp = resp / max(resp);
subplot 313, plot(resp),title('Response to the whole set of stimuli')

%%  hrf convolved multiple stimuli
res = 100;
% this is the number of samples pre second - the temporal resolution
% when creating the time course. We'll downsample it to match what the
% scanner can really do later.
duration = 200; % seconds
% create a canonical BOLD response - impulse response function
hrf = spm_hrf( 1/res );
% create two waveforms representing the activity of two mental processes
stim1 =zeros(duration*res,1);
stim1( 5*res:13*res:end ) =1;
stim2 =zeros(duration*res,1);
for n=8:10:duration-15
 stim2( n*res : (n+4)*res ) =1;
end
figure(5)
subplot(311)
plot(stim1); hold on

axis([ 0 duration*res -1 2])
subplot(312)
plot(stim2); hold on
axis([ 0 duration*res -1 2])
title('Two types of events')
%Next we convolve the input and the HRF
% first note what the convolution operation does
resp1 = conv(stim1, hrf);
resp2 = conv(stim2, hrf);
% Note the size s of the inputs and the outputs to the function:
whos resp stim hrf
% remove the data at the end of the response function, so that it matches the
% size of the data, and normalize it
resp1 = resp1(1:length(stim1));
resp1 = resp1 / max(resp1);
resp2 = resp2(1:length(stim2));
resp2 = resp2 / max(resp2);
subplot 311,
plot(resp1,'r'),hold off
title('Event type 1')
subplot 312,
plot(resp2,'r'),hold off
title('Event Type 2')

% More realistically, the two processes will rarely have the same amplitude, but the observed signal will be a
% combination of the two processes with different amplitudes, something like this:
% Now let's put these two things together into a single time course.
% Each of the responses will have different amplitude - let's call that
% beta
beta2 = 0.5;
beta1 = 1;
beta0 = 10; % this will be the baseline signal
% This is the simple way to construct a model:
% you add the different components of the signal
% weighted by some coefficient:
y = beta0 + beta1 * resp1 + beta2*resp2;
% You?ll note that there is an offset added to the signal (beta0). This is because there is always some
% baseline signal present. This is just the image intensity of the MRI image and is not related to brain
% activity.
subplot(313)
plot(y); title ('The whole signal ')

