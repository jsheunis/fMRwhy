% fmrwhy_script_createNewTmapMontages

% Main input: BIDS root folder
bids_dir = '/Volumes/TSM/NEUFEPME_data_BIDS';

all_reports_dir = '/Users/jheunis/Desktop/multiecho_reports2/';

% Loop through subjects
subs = {'002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};
%subs = {'001'};

for s = 1:numel(subs)
    sub = subs{s};
    reports_dir = fmrwhy_neufep_generateMEsubReport(bids_dir, sub);
    copyfile(reports_dir, [all_reports_dir 'sub-' sub])
end