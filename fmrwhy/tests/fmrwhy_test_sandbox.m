% p = mfilename('fullpath')
% IND = strfind(p,'fMRwhy')
% fmrwhy_dir = p(1:IND+5)
%
%
% src = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-qc/sub-001/sub-001_task-emotion_run-1_desc-QCreport_2020430225716.html';
% dest = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-qc/sub-001/report/koekier.html';
% copyfile(src, dest)

bids_dir = '/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS';
sub = '001';
fmrwhy_neufep_generateSubReport(bids_dir, sub);
