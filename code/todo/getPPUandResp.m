function [x, ppu, resp, outParams, waveforms] = getPPUandResp(logfile, TR, Ndyn)

% TR in ms

% fn2 = [template_dir filesep subj filesep 'scanphyslog/SCANPHYSLOG20170816100312.log'];

% From welcheb function:
[outParams, waveforms] = loadSCANPHYSLOG(logfile)
f_sampling = 500; % Sampling rate of Philips physiology sensors (if hardwired) 

start_ind = find(waveforms.mark_dec == 256);
end_ind = find(waveforms.mark_dec == 32);
end_ind = end_ind(end);
% time2 = (end_ind - start_ind)/500;
t_sampling = 1/f_sampling;
ind_end = end_ind;
length = (Ndyn*(TR/1000)*f_sampling);
ind_start = ind_end - length + 1;
resp = waveforms.resp(ind_start:ind_end);
ppu = waveforms.ppu(ind_start:ind_end);
x = (1:length)/1000; 
end