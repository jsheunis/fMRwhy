function output = fmrwhy_qc_createSubReport(bids_dir, sub, ses, task, run)

% Load/create required defaults
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
template_task = 'motor'; % changed for fingertapping experiment. TODO: change back. and update functioning.
template_run = '1';
template_echo = '2';
stre_txt = ['sub-' sub '_task-' task '_run-' run '_echo-' echo];
str_txt = ['sub-' sub '_task-' task '_run-' run];

% Directory and content setup
deriv_dir = fullfile(bids_dir, 'derivatives');
preproc_dir = fullfile(deriv_dir, 'fmrwhy-preproc');
qc_dir = fullfile(deriv_dir, 'fmrwhy-qc');
sub_dir_preproc = fullfile(preproc_dir, ['sub-' sub]);
sub_dir_qc = fullfile(qc_dir, ['sub-' sub]);
func_dir_qc = fullfile(sub_dir_qc, 'func');
anat_dir_qc = fullfile(sub_dir_qc, 'anat');


% Get sample HTML file
fid = fopen('in.txt');
lines = textscan(fid,'%s','delimiter','\n');
fclose(fid);
lines = lines{1};
%replace tags with actual values
relevant =  find(~cellfun(@isempty,strfind(s{1},'E_FROM')));
for i = 1:length(relevant)
    lines{relevant(i)-14} = '<modified line>';
end
% Create subject report
fid = fopen('out.txt','w');
for i = 1:length(lines)
    fprintf(fid,'%s\n',lines{i});
end
fclose(fid)






















