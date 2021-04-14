function noMoxTest_fmrwhy_qc_calculateFD

    %% define the context of the test

    % define the path to test data
    % relative to the path of the function we are running
    test_data_dir = fullfile(fileparts(mfilename('fullpath')), 'test_data');

    movement_parameter_file = fullfile(test_data_dir, 'rp_sub-01_task-auditory_bold.txt');

    movement_parameter = spm_load(movement_parameter_file);

    radius = 60;

    FD_threshold = 0;

    %% run the function  we are tesring

    FD_measures = fmrwhy_qc_calculateFD(movement_parameter, radius, FD_threshold);

    %% define the data we want to compare our ouput to

    expected_data_file = fullfile(test_data_dir, 'expected_FD_thres_0.mat');

    % when creating test you can use this to create the expected output
    %
    %     expected_FD = FD_measures;
    %
    %     save(expected_data_file, 'expected_FD')

    load(expected_data_file);

    %% run the test comparing that the output we get is equal to what we expect

    % we can do it on the whole structure
    assert(isequal(FD_measures, expected_FD));

    % or on some fields
    assert(isequal(FD_measures.FD_sum, expected_FD.FD_sum));

    % if everything is equal then we are all good and no error is thrown

end
