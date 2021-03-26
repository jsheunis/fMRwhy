% This script will download the dataset from the FIL for the block design SPM tutorial
% and will run the fMRwhy on it.

working_directory = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(working_directory, '..', '..', 'fmrwhy')));
addpath(genpath(fullfile(working_directory, '..', '..', 'lib')));

% download_data = true;
dowload_MoAE_ds(download_data);

fmrwhy_workflow_qc(fullfile(working_directory, 'MoAE_settings.m'));



%%
function dowload_MoAE_ds(downloadData)

  if downloadData

    % URL of the data set to download
    URL = 'http://www.fil.ion.ucl.ac.uk/spm/download/data/MoAEpilot/MoAEpilot.bids.zip';

    working_directory = fileparts(mfilename('fullpath'));

    % clean previous runs
    if exist(fullfile(working_directory, 'data'), 'dir')
      rmdir(fullfile(working_directory, 'data'), 's');
    end

    spm_mkdir(fullfile(working_directory, 'data'));

    %% Get data
    fprintf('%-10s:', 'Downloading dataset...');
    urlwrite(URL, 'MoAEpilot.zip');
    fprintf(1, ' Done\n\n');

    fprintf('%-10s:', 'Unzipping dataset...');
    unzip('MoAEpilot.zip');
    movefile('MoAEpilot', fullfile(working_directory, 'data'));
    fprintf(1, ' Done\n\n');

  end

end



