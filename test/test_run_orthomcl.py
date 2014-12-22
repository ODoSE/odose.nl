from Bio import SeqIO
import logging
import os
import tempfile
import unittest

import run_orthomcl
from shared import resource_filename
from versions import ORTHOMCL_DIR


class Test(unittest.TestCase):

    def setUp(self):
        self.longMessage = True
        logging.root.setLevel(logging.DEBUG)

    def test_parse_args(self):
        args = run_orthomcl._parse_args(['proteins.zip', '-e', '4', 'target-poor.fasta', 'target-groups.tsv'])
        logging.debug(args)
        self.assertEqual(4, args.evalue)
        self.assertEqual(30, args.poorlength)

    @unittest.skipUnless(os.path.isdir(ORTHOMCL_DIR), 'We need OrthoMCL')
    def test_run_orthomcl(self):
        '''
        Run run_orthomcl.run_orthomcl on two genomes and verify the poor proteins and identified groups
        '''
        # Setup
        proteome_files = [resource_filename(__name__, 'data/run_orthomcl/' + acc + '.1.faa') for acc in ['13305', '17745']]
        poor_protein_length = 30
        evalue_exponent = -5
        target_poor_proteins_file = tempfile.mkstemp(suffix='.txt', prefix='poor_proteins_')[1]
        target_groups_file = tempfile.mkstemp(suffix='.txt', prefix='groups_')[1]

        try:
            # Exercise
            run_orthomcl.run_orthomcl(proteome_files, poor_protein_length, evalue_exponent, target_poor_proteins_file, target_groups_file)

            # Verify
            # Poor proteins that must be skipped
            poor_seqr = [rec for rec in SeqIO.parse(target_poor_proteins_file, 'fasta')]
            self.assertEqual('13305_1|YP_667942.1', poor_seqr[0].id)
            self.assertEqual('13305_1|YP_668016.1', poor_seqr[1].id)

            # Groups that must be identified
            groups = sorted([['17745_1|YP_001451418.1', '17745_1|YP_001451593.1', '13305_1|YP_667991.1'],
                             ['17745_1|YP_001451427.1', '17745_1|YP_001451583.1', '13305_1|YP_668281.1'],
                             ['17745_1|YP_001451428.1', '17745_1|YP_001451499.1', '13305_1|YP_671707.1'],
                             ['17745_1|YP_001451436.1', '17745_1|YP_001451438.1', '13305_1|YP_668218.1'],
                             ['13305_1|YP_670193.1', '13305_1|YP_671256.1', '17745_1|YP_001451410.1'],
                             ['13305_1|YP_670987.1', '13305_1|YP_671559.1', '13305_1|YP_672523.1'],
                             ['13305_1|YP_670989.1', '13305_1|YP_671558.1', '13305_1|YP_672521.1'],
                             ['13305_1|YP_671185.1', '13305_1|YP_672343.1', '13305_1|YP_671555.1'],
                             ['13305_1|YP_671533.1', '13305_1|YP_672559.1', '13305_1|YP_671959.1'],
                             ['13305_1|YP_671597.1', '13305_1|YP_671599.1', '13305_1|YP_671600.1'],
                             ['17745_1|YP_001451452.1', '17745_1|YP_001451477.1'],
                             ['17745_1|YP_001451451.1', '17745_1|YP_001451476.1'],
                             ['17745_1|YP_001451453.1', '17745_1|YP_001451478.1'],
                             ['17745_1|YP_001451458.1', '17745_1|YP_001451482.1'],
                             ['17745_1|YP_001451460.1', '17745_1|YP_001451468.1'],
                             ['13305_1|YP_672143.1', '13305_1|YP_672524.1'],
                             ['13305_1|YP_672210.1', '13305_1|YP_672227.1'],
                             ['13305_1|YP_672222.1', '13305_1|YP_672446.1'],
                             ['13305_1|YP_672237.1', '13305_1|YP_672238.1'],
                             ['13305_1|YP_672511.1', '13305_1|YP_672516.1'],
                             ['13305_1|YP_671068.1', '17745_1|YP_001451503.1'],
                             ['13305_1|YP_671069.1', '17745_1|YP_001451504.1'],
                             ['13305_1|YP_671072.1', '17745_1|YP_001451505.1'],
                             ['13305_1|YP_671147.1', '17745_1|YP_001451543.1'],
                             ['13305_1|YP_671255.1', '17745_1|YP_001451409.1'],
                             ['13305_1|YP_671502.1', '17745_1|YP_001451590.1'],
                             ['13305_1|YP_671676.1', '17745_1|YP_001451378.1'],
                             ['13305_1|YP_671710.1', '17745_1|YP_001451500.1'],
                             ['13305_1|YP_672136.1', '17745_1|YP_001451472.1'],
                             ['13305_1|YP_672448.1', '17745_1|YP_001451386.1']])
            with open(target_groups_file) as reader:
                contents = reader.read()
                actual = sorted([line.split('\t') for line in contents.split('\n') if line])
            self.assertEqual(groups, actual)

        finally:
            os.remove(target_groups_file)
            os.remove(target_poor_proteins_file)

    @unittest.skipUnless(os.path.isdir(ORTHOMCL_DIR), 'We need OrthoMCL')
    def test_main(self):
        proteins = resource_filename(__name__, 'data/run_orthomcl/proteins.zip')
        poor = tempfile.mkstemp(suffix='.faa', prefix='poor_')[1]
        groups = tempfile.mkstemp(suffix='.tsv', prefix='groups_')[1]
        try:
            run_orthomcl.main([proteins, poor, groups])
        finally:
            os.remove(poor)
            os.remove(groups)
