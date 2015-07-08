'''
Created on Jun 8, 2014

@author: tim
'''
import datetime
import logging
import os
import sys
import tempfile
import unittest

import load_prokaryotes


class Test(unittest.TestCase):

    # Assertions about particular record
    # Escherichia coli E24377A    331111    PRJNA13960    13960    Proteobacteria    Gammaproteobacteria    5.24929    50.5414    NC_009801.1    CP000800.1    NC_009786.1,NC_009789.1,NC_009788.1,NC_009787.1,NC_009790.1,NC_009791.1    CP000795.1,CP000798.1,CP000797.1,CP000796.1,CP000799.1,CP000801.1    -    7    5258    4991    2007/09/11    2014/01/31    Gapless Chromosome    TIGR    SAMN02604038    GCA_000017745.1    -    Escherichia_coli/GCF_000017745    18676672
    ref = {'Assembly Accession': '17745.1',
           'BioProject Accession': 'PRJNA13960',
           'BioProject ID': '13960',
           'BioSample Accession': 'SAMN02604038',
           'Center': 'TIGR',
           'Chromosomes/INSDC': ['CP000800.1'],
           'Chromosomes/RefSeq': ['NC_009801.1'],
           'FTP Path': 'Escherichia_coli/GCF_000017745',
           'GC%': '50.5414',
           'Genes': '5258',
           'Group': 'Proteobacteria',
           'Modify Date': datetime.datetime(2014, 1, 31, 0, 0),
           'Organism/Name': 'Escherichia coli O139:H28 str. E24377A',
           'Plasmids/INSDC': ['CP000795.1',
                              'CP000798.1',
                              'CP000797.1',
                              'CP000796.1',
                              'CP000799.1',
                              'CP000801.1'],
           'Plasmids/RefSeq': ['NC_009786.1',
                               'NC_009789.1',
                               'NC_009788.1',
                               'NC_009787.1',
                               'NC_009790.1',
                               'NC_009791.1'],
           'Proteins': '4991',
           'Pubmed ID': '18676672',
           'Reference': '-',
           'Release Date': datetime.datetime(2007, 9, 11, 0, 0),
           'Scaffolds': '7',
           'Size (Mb)': '5.24929',
           'Status': 'Gapless Chromosome',
           'SubGroup': 'Gammaproteobacteria',
           'TaxID': '331111',
           'WGS': '-'}

    def setUp(self):
        self.longMessage = True
        logging.root.setLevel(logging.DEBUG)

    def test_parse_genomes_table(self):
        # Parse all genomes
        genomes = load_prokaryotes._parse_genomes_table()
        self.assertLess(2500, len(genomes), 'Enough genomes should be available')

        # Determine the unique columns
        keys = genomes[0].keys()
        uniq = [key for key in keys if len(set(str(gnm[key]) for gnm in genomes)) == len([gnm[key] for gnm in genomes])]
        logging.debug('Unique columns: %s', uniq)

        # Ensure the Assembly accessions are unique
        self.assertIn('Assembly Accession', uniq)

        # Report how many records RefSeq / GenBank identifiers
        self.assertTrue(all(gnm['Chromosomes/RefSeq'] for gnm in genomes))

        # Ensure there are no duplicates in Chromosomes/RefSeq
        for gnm in genomes:
            self.assertEqual(len(gnm['Chromosomes/RefSeq']), len(set(gnm['Chromosomes/RefSeq'])))

        # Extract a particular record
        for gnm in genomes:
            if gnm['BioProject ID'] == '13960':
                print gnm
                break
        else:
            self.fail('Escherichia coli E24377A not found!')

        for prop in ['Organism/Name', 'FTP Path', 'BioProject ID', 'TaxID', 'Assembly Accession']:
            self.assertEqual(self.ref[prop], gnm[prop])

    def test_main(self):
        # Setup arguments
        target = tempfile.mktemp()[1]
        try:
            sys.argv = ['', target]

            # Write to argument file
            load_prokaryotes.main()

            # Assert contents
            with open(target) as reader:
                contents = reader.read()
            self.assertIn('439255.1 - Proteobacteria > Gammaproteobacteria > Salmonella > Salmonella bongori N268-08', contents)
        finally:
            os.remove(target)
