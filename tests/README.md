# Tests for fMRwhy

We use a series of unit and integration tests to make sure the `fMRwhy` code behaves as
expected and to also help in development.

If you are not sure what unit and integration tests are, check out the excellent
chapter about these topics in the
[Turing way](https://the-turing-way.netlify.app/reproducible-research/testing.html) documentation.

## How to run the tests

### 0. Example test

`fMRwhy` uses [MOxUnit for matlab and octave](https://github.com/MOxUnit/MOxUnit) to run unit and integration tests,
but a separate example test script is also provided.
This test script, `noMoxTest_fmrwhy_qc_calculateFD.m` (located in the `fMRwhy/tests` directory),
can be run directly from the MATLAB / Octave command window without having MOxUnit installed.
It serves as a minimal example of how a unit test can be constructed for `fMRwhy` code.

Test data (located in `fMRwhy/tests/test_data`) is also provided along with the example script.
This data is also used by the MOxUnit test suite.


### 1. Install MOxUnit

You need to install
[MOxUnit for matlab and octave](https://github.com/MOxUnit/MOxUnit) to run the
test suite.

Note, the install procedure will require you to have
[git](https://git-scm.com/downloads) installed on your computer. If you don't,
you can always download the MOxUnit code with this
[link](https://github.com/MOxUnit/MOxUnit/archive/master.zip).

Run the following from your terminal after navigating to the directory where you want to install
MOxUnit. The `make install` command will find MATLAB / Octave on your system and
make sure it plays nice with MOxUnit.

```bash
# get the code for MOxUnit with git
git clone https://github.com/MOxUnit/MOxUnit.git
# enter the newly created folder and set up MoxUnit
cd MOxUnit
make install
```

If you want to check the code coverage on your computer, you can also install
[MOcov for matlab and octave](https://github.com/MOcov/MOcov). Note that this package is
also part of the continuous integration setup of the `bids-matlab` package,
which in turn is a dependency of `fMRwhy`. So you might already have MOcov installed on your system,
depending on how you installed `bids-matlab`.

### 3. Run the tests

From the root folder of the `fMRwhy` directory, you can run the test suite with one
the following commands.

```bash
moxunit_runtests tests
# Or if you want more feedback
moxunit_runtests tests -verbose
```

MOxUnit will automatically detect all tests that need to be run (currently contained in `test_fmrwhy_qc_calculateFD.m`).
Test results and other outputs will be displayed in the MATLAB / Octave command window.


## Adding more tests

You can use the following function template to write more tests.

```matlab
function test_suite = test_functionToTest()
    % This top function is necessary for mox unit to run tests.
    % DO NOT CHANGE IT except to adapt the name of the function.
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions = localfunctions(); %#ok<*NASGU>
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
end

function test_function_to_test_basic()

    %% set up

    %% data to test against

    %% test
    % assertTrue( );
    % assertFalse( );
    % assertEqual( );

end

function test_function_to_test_other_usecase()

    %% set up

    %% data to test against

    %% test
    % assertTrue( );
    % assertFalse( );
    % assertEqual( );

end

```
