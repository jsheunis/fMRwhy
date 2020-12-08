% download open dataset from Marcus Munafo 
% https://data.bris.ac.uk/datasets/112g2vkxomjoo1l26vjmvnlexj/2016.08.14_AnxietyPaper_Data%20Sheet.csv

% % mydir = 'insert path here';
% mydir = '/Users/jheunis/Documents/MATLAB';
% myfile = '2016.08.14_AnxietyPaper_Data Sheet.csv';
% 
% f = fopen(fullfile(mydir,myfile));
% 
% dall = textscan(f,'','HeaderLines',1,'delimiter',',');
% 
% % just get the Anger, Disgust, Fear, and Happy conditions
% d = [dall{28} dall{29} dall{30} dall{31}];
% 
% fclose(f);

%% get nice colours from color brewer
% (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)

[cb] = cbrewer('qual','Set3',10,'pchip');

%% plot

figure;
raincloud_plot(FD_measures.FD, cb(5,:))
view([90 90])
set(gca, 'Xdir', 'reverse');

% subplot(4,1,1), raincloud_plot(d(:,1), cb(5,:));
% subplot(4,1,2), raincloud_plot(d(:,2), cb(7,:));
% subplot(4,1,3), raincloud_plot(d(:,3), cb(6,:));
% subplot(4,1,4), raincloud_plot(d(:,4), cb(4,:));