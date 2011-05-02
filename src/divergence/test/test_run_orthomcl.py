#!/usr/bin/env python
"""Test module for divergence.run_orthomcl"""

from divergence import resource_filename
from divergence.run_orthomcl import run_orthomcl
from divergence.test import fail_if_diff
import os
import unittest

class Test(unittest.TestCase):
    """Test class"""

    def test_minor_unrelated_set(self):
        """Run orthomcl for two sequences to pollute files or database."""
        aa_files = []
        aa_files.append(resource_filename(__name__, 'data/orthomcl/57915.faa'))
        aa_files.append(resource_filename(__name__, 'data/orthomcl/58017.faa'))

        #Run orthomcl
        groups_file = run_orthomcl(aa_files)
        self.failUnless(os.path.isfile(groups_file))
        self.failUnless(0 < os.path.getsize(groups_file))

        #Compare actual output file with expected output file
        expected_groups_file = resource_filename(__name__, 'data/orthomcl/groups-minor-set.txt')
        fail_if_diff(self, expected_groups_file, groups_file)

    def test_orthomcl_ecoli_salmo(self):
        """Run orthomcl for first two refseq E coli vs first two refseq Salmonella."""
        aa_files = []
        aa_files.append(resource_filename(__name__, 'data/orthomcl/58191.faa'))
        aa_files.append(resource_filename(__name__, 'data/orthomcl/58531.faa'))
        aa_files.append(resource_filename(__name__, 'data/orthomcl/59245.faa'))
        aa_files.append(resource_filename(__name__, 'data/orthomcl/59431.faa'))

        #Run orthomcl
        groups_file = run_orthomcl(aa_files)
        self.failUnless(os.path.isfile(groups_file))
        self.failUnless(0 < os.path.getsize(groups_file))

        #Compare actual output file with expected output file
        expected_groups_file = resource_filename(__name__, 'data/orthomcl/groups.txt')
        fail_if_diff(self, expected_groups_file, groups_file)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
