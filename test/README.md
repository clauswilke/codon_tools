# Writing unit tests

Unit tests are auto-discovered based on the module name. So, unit tests for `codon_analyzer.py` should go in `codon_analyzer_test.py`. The test runner will always attempt to match a `*_test.py` filename pattern.

All unit tests should extend `unittest.TestCase`. The `unittest` package provides helper functions for setting up code and cleaning up code after tests are run. They are not required but can be useful.

Within a `unittest.TestCase` child class, the test runner will run any method with the prefix `test_`. Within each of these test methods, you can call many different assertion methods.

Run: `setup.py test` to see test output.

