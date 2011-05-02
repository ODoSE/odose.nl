#!/usr/bin/env python
"""Module to test translation between DNA and protein level."""

from divergence import resource_filename, create_directory
from divergence.translate import _extract_gene_and_protein, _build_protein_to_cog_mapping, translate_genomes, main
from divergence.test import ECOLI_GENOMES, SALMONELLA_GENOMES, fail_if_diff
import os.path
import unittest

class Test(unittest.TestCase):
    """Test class"""

    def test_protein_to_cog_mapping(self):
        """Test mapping protein to COG."""
        ptt_file = resource_filename(__name__, 'data/translate/NC_012947.ptt')
        mapping = _build_protein_to_cog_mapping(ptt_file)
        #for key in sorted(mapping.keys()):
        #    print(key + ': ' + str(mapping[key]))
        self.failUnlessEqual(4228, len(mapping))

    def test_translate_genbank_to_prot(self):
        """Test translation from DNA to protein."""
        genbank_file = resource_filename(__name__, 'data/translate/NC_012947.gbk')
        ptt_file = resource_filename(__name__, 'data/translate/NC_012947.ptt')

        test_dir = create_directory('translations')
        dna_file, aa_file = _extract_gene_and_protein(test_dir, 'testid', genbank_file, ptt_file)

        self.assertTrue(os.path.exists(dna_file))
        self.assertTrue(os.path.exists(aa_file))

        exp_dna_file = resource_filename(__name__, 'data/translate/NC_012947.ffn')
        exp_aa_file = resource_filename(__name__, 'data/translate/NC_012947.faa')

        fail_if_diff(self, exp_dna_file, dna_file)
        fail_if_diff(self, exp_aa_file, aa_file)

    def test_with_known_reference(self):
        """Test translation of some known files."""
        genomes = ECOLI_GENOMES[:2] + SALMONELLA_GENOMES[:2]
        dna_files, protein_files = translate_genomes(genomes)
        for dna_file in dna_files:
            self.failUnless(0 < os.path.getsize(dna_file))
            refseq_id = os.path.splitext(os.path.split(dna_file)[1])[0]
            fail_if_diff(self, resource_filename(__name__, 'data/orthomcl/{0}.ffn'.format(refseq_id)), dna_file)
        for protein_file in protein_files:
            self.failUnless(0 < os.path.getsize(protein_file))
            refseq_id = os.path.splitext(os.path.split(protein_file)[1])[0]
            fail_if_diff(self, resource_filename(__name__, 'data/orthomcl/{0}.faa'.format(refseq_id)), protein_file)

    def test_main(self):
        """Test running main method of translate.py by building a """
        #Determine in- and output filenames
        ids_file = resource_filename(__name__, 'data/translate/genome_ids.txt')
        target_dna = '/tmp/dna.zip'
        target_protein = '/tmp/protein.zip'

        #Run main with below arguments
        args = ['--genomes', ids_file, '--dna-zip', target_dna, '--protein-zip', target_protein]
        main(args)

        #Assert output file paths now exist
        self.failUnless(0 < os.path.getsize(target_dna))
        self.failUnless(0 < os.path.getsize(target_protein))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
