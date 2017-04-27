import unittest
from codon_tools.codon_analyzer import CodonAnalyzer

'''
Writing unit tests

Unit tests are auto-discovered based on the module name. So, unit tests for 'codon_analyzer.py' should go in 'codon_analyzer_test.py'. The test runner will always attempt to match a '*_test.py' filename pattern.

All unit tests should extend unittest.TestCase. The unittest package provides
helper functions for setting up code and cleaning up code after tests are run. They are not required but can be useful.

Within a unittest.TestCase child class, the test runner will run any method with the prefix 'test_'. Within each of these test methods, you can call many different assertion methods. I show a few examples in the code below.

Run: `setup.py test` to see test output.
'''

class TestCodonAnalyzer(unittest.TestCase):
    def setUp(self):
        # Do any set up that you need here
        self.my_analyzer = CodonAnalyzer()
    
    def test_count_stop_muts(self):
        seq1 = "ATG"
        seq2 = "CGT"
        seq3 = "AGA"

        # Here are some tests
        self.assertEqual(self.my_analyzer.count_stop_muts(seq1), 0)
        self.assertEqual(self.my_analyzer.count_stop_muts(seq2), 0)
        self.assertEqual(self.my_analyzer.count_stop_muts(seq3), 1)

        # Here are some meaningless test examples that demonstrate assertions
        self.assertTrue(1 == 1)
        self.assertFalse(1 == 0)
        self.assertEqual("hello", "hello")
        self.assertIsNone(None)
        # You can also test that exceptions are raised
        with self.assertRaises(AssertionError):
            # Code that should raise an error here
            self.assertEqual(self.my_analyzer.count_stop_muts(seq1 + seq2 + seq3), 3)
    
    def tearDown(self):
        # Do any clean up here
        pass

