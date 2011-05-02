#!/usr/bin/env python
"""Test module for entire analysis"""

from divergence import resource_filename
from divergence.clean_orthologs import trim_and_concat_sicos
from divergence.run_orthomcl import run_orthomcl
from divergence.run_paml import run_paml
from divergence.select_taxa import select_taxa
from divergence.test import fail_if_diff
from divergence.translate import translate_genomes
from divergence.test import ECOLI_IDS, SALMO_IDS
import os.path
import unittest

class Test(unittest.TestCase):
    """Test class"""

    def test_all_ecoli_salmo(self):
        """Run entire analysis for selected taxa."""
        #Let user select taxa
        genomes_a, genomes_b = select_taxa(ECOLI_IDS, SALMO_IDS)

        #Subset genomes to make the test run faster
        genomes_a, genomes_b = genomes_a[-2:], genomes_b[:2]
        genomes = genomes_a + genomes_b

        #Download & translate
        dna_files, protein_files = translate_genomes(genomes)

        #Run orthomcl
        groups_file = run_orthomcl(protein_files)

        #Run post orthomcl cleanup steps
        trimmed_sico_files, concatemers, stats_file = trim_and_concat_sicos(genomes, dna_files, groups_file)

        for sico_file in trimmed_sico_files:
            self.failUnless(0 < os.path.getsize(sico_file), 'Each trimmed alignment file should exists with content')
        for genome_file in concatemers:
            self.failUnless(0 < os.path.getsize(genome_file), 'Each genome concatemer should exists with content')
        self.failUnless(0 < os.path.getsize(stats_file), 'Stats file should exist with content')

        #Diff actual outputs with expected outputs
        expected_58191 = resource_filename(__name__, 'data/test_all/58191.trimmed.concatemer.fasta')
        fail_if_diff(self, expected_58191, sorted(concatemers)[0])

        expected_stats = resource_filename(__name__, 'data/test_all/stats.txt')
        fail_if_diff(self, expected_stats, stats_file)

        #Run PAML to calculate dN/dS per SICO
        run_paml(genomes_a, genomes_b, trimmed_sico_files)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
