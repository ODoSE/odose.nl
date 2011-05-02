#!/usr/bin/env python
"""Test module for divergence.select_taxa"""

from divergence import resource_filename
from divergence.select_taxa import _parse_genomes_table, _print_genomes_tree, select_taxa, download_genome_files
from divergence.test import ECOLI_IDS, SALMO_IDS
import os
import unittest

class Test(unittest.TestCase):
    """Test class"""
    def test_parse_genomes_table(self):
        """Parse table of complete microbial genomes and print returned list and dictionaries of genomes."""
        with open(resource_filename(__name__, 'data/genomes_table.tsv')) as table_file:
            content = table_file.read().decode('utf-8')
        refseq_genomes = _parse_genomes_table(content)
        self.failUnless('RefSeq project ID' in refseq_genomes[0])
        self.failUnless('List of GenBank accessions' in refseq_genomes[0])
        _print_genomes_tree(refseq_genomes)

        #Parse again, but now do not filter out genomes without refseq ids
        genbank_genomes = _parse_genomes_table(content, require_refseq = False)
        self.assertTrue(len(refseq_genomes) < len(genbank_genomes), 'More genomes should be found when req is dropped.')

    def test_download_selected_taxa(self):
        """Select four taxa from list of genomes and assert they match expected values."""
        clade_a_ids = ECOLI_IDS[:2]
        clade_b_ids = SALMO_IDS[:2]
        clade_a, clade_b = select_taxa(clade_a_ids, clade_b_ids)

        self.failUnlessEqual(set(clade_a_ids), set([genome['RefSeq project ID'] for genome in clade_a]), clade_a_ids)
        self.failUnlessEqual(set(clade_b_ids), set([genome['RefSeq project ID'] for genome in clade_b]), clade_b_ids)

        #Assert files are translated correctly
        for genome in clade_a + clade_b:
            for gbk_file, ptt_file in download_genome_files(genome):
                self.failUnless(os.path.isfile(gbk_file) and 0 < os.path.getsize(gbk_file), gbk_file)
                self.failUnless(os.path.isfile(ptt_file) and 0 < os.path.getsize(ptt_file), ptt_file)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
