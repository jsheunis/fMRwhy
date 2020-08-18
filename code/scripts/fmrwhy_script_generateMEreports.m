% fmrwhy_script_createNewTmapMontages

% Main input: BIDS root folder
bids_dir = '/Users/jheunis/Desktop/NEUFEPME_data_BIDS';
% Loop through subjects
subs = {'001', '002', '003', '004', '005', '006', '007', '010', '011', '012', '013', '015', '016', '017', '018', '019', '020', '021', '022', '023', '024', '025', '026', '027', '029', '030', '031', '032'};

for s = 1:numel(subs)
    sub = subs{s};
    fmrwhy_neufep_generateMEsubReport(bids_dir, sub)
end