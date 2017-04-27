import unittest
from codon_tools.codon_analyzer import CodonAnalyzer


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

