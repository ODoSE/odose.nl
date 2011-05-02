#!/usr/bin/env python
"""Test module for divergence.run_paml"""

from divergence import resource_filename
from divergence.run_paml import run_paml
from divergence.test import ECOLI_GENOMES, SALMONELLA_GENOMES
import unittest

class Test(unittest.TestCase):
    """Test class"""

    def test_run_paml(self):
        """Run PAML for 10 SICO files."""
        #Select genomes matching SICO test files
        genomes_a, genomes_b = ECOLI_GENOMES[:2], SALMONELLA_GENOMES[:2]

        #Retrieve test sico files
        sico_files = []
        sico_files.append(resource_filename(__name__, 'data/paml/sico_000000.nt_ali.fasta'))
        sico_files.append(resource_filename(__name__, 'data/paml/sico_000001.nt_ali.fasta'))
        sico_files.append(resource_filename(__name__, 'data/paml/sico_000002.nt_ali.fasta'))
        sico_files.append(resource_filename(__name__, 'data/paml/sico_000003.nt_ali.fasta'))
        sico_files.append(resource_filename(__name__, 'data/paml/sico_000004.nt_ali.fasta'))

        #Run paml
        run_paml(genomes_a, genomes_b, sico_files)
        self.failUnless(True)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
