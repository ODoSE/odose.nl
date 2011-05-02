#!/usr/bin/env python
"""Test module for divergence.reciprocal_blast"""

from divergence import resource_filename
from divergence.reciprocal_blast import reciprocal_blast
from divergence.test import fail_if_diff
import unittest

class Test(unittest.TestCase):
    """Test class"""

    def test_reciprocal_blast(self):
        """Blast a small subset of the kamchatka & california sulfolobus proteomes against eachother."""
        tx1_and_tx2 = resource_filename(__name__, 'data/blast/tx1_and_tx2.fasta')
        tx1 = resource_filename(__name__, 'data/blast/tx1.fasta')
        tx2 = resource_filename(__name__, 'data/blast/tx2.fasta')
        allvsallhits = reciprocal_blast(tx1_and_tx2, [tx1, tx2])
        expected_hits = resource_filename(__name__, 'data/blast/expected_all-vs-all.tsv')
        fail_if_diff(self, expected_hits, allvsallhits)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
