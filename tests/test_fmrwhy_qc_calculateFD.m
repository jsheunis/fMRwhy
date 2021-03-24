function test_suite = test_fmrwhy_qc_calculateFD %#ok<*STOUT>
  try % assignment of 'localfunctions' is necessary in Matlab >= 2016
    test_functions = localfunctions(); %#ok<*NASGU>
  catch % no problem; early Matlab versions can use initTestSuite fine
  end
  initTestSuite;
end

function test_fmrwhy_qc_calculateFD_no_threshold()
    
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
    
    load(expected_data_file)
    
    %% use assertEqual function from MoxUnit as it is better suited than a simple assert
   
    assertEqual(FD_measures, expected_FD);
    
    % if everything is equal then we are all good and no error is thrown

end

function test_fmrwhy_qc_calculateFD_with_threshold()
    
     %% define the context of the test
    
    % define the path to test data
    % relative to the path of the function we are running
    
    % TODO this few lines should probably be refactored into a little helper
    % function that just loads the data and does not distract us from the
    % content of the actual test
    
    test_data_dir = fullfile(fileparts(mfilename('fullpath')), 'test_data');
    
    movement_parameter_file = fullfile(test_data_dir, 'rp_sub-01_task-auditory_bold.txt');
    
    movement_parameter = spm_load(movement_parameter_file);
    
    radius = 60;
    
    FD_threshold = 2;
    
    
    %% run the function  we are tesring
    
    FD_measures = fmrwhy_qc_calculateFD(movement_parameter, radius, FD_threshold);
    
    
    %% define the data we want to compare our ouput to
    
    expected_data_file = fullfile(test_data_dir, 'expected_FD_thres_0.mat');
    
    % when creating test you can use this to create the expected output
    %
    %     expected_FD = FD_measures;
    %     
    %     save(expected_data_file, 'expected_FD')
    
    load(expected_data_file)
    
    %% use assertEqual function from MoxUnit as it is better suited than a simple assert

    % ----------------------- % !!!!!!!!!!!!!!!!!!!!! % ----------------------- %
    
    %  THIS WILL FAIL BECAUSE WE DON'T HAVE TEST DATA FOR --> FD_threshold = 2;
    
    assertEqual(FD_measures, expected_FD);
    
    % ----------------------- % !!!!!!!!!!!!!!!!!!!!! % ----------------------- %

end