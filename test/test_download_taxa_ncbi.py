'''
Created on Jun 8, 2014

@author: tim
'''
import datetime
import os
import tempfile
import unittest

import download_taxa_ncbi


class Test(unittest.TestCase):

    def test_download_genome_files(self):
        # Example genome record extracted once from select_taxa
        # Maps to: ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_E24377A_uid58395/
        gnm = {'Assembly Accession': '17745.1',
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
               'Organism/Name': 'Escherichia coli E24377A',
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
        log= tempfile.mkstemp()[1]
        try:
            download_taxa_ncbi.download_genome_files(gnm, download_log=log)
            with open(log) as reader:
                logcontent = reader.read()
        finally:
            os.remove(log)
        self.assertNotIn('Genome skipped because of missing files', logcontent)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
