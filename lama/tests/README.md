# LAMA tests

There are a number of tests included in LAMA, located at lama/tests
The lama suite of tests is currently under development and most of the tests avaialble jsut run a part of the pipeline
and checks whether it runs to completion without exceptions being raised.

### Test data
To run the test, you must download the test data. These data consist of small volumes for registration and config files 
that are designed to run rapidly to make sure the components of thepipeline are working correctly. 
The tests are not made to test the accuracy of the registration.

```bash
lama_get_test_data
```
Or download it form here
[[test data|http://images.mousephenotype.org/lama/test_data.zip]]
Uzip the test_data folder and place in  /lama/tests

### Run the tests
Test are run using [[PyTest|https://docs.pytest.org]]

```bash
lama_phenotype_detection/tests/run_tests.sh
```


### TODO
* Make more tests
* Make tests that query the output data to ensure everything has been procesed correctly
* Make the permutation tests output to main stats folder
