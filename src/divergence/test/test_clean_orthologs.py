#!/usr/bin/env python
"""Test module for divergence.clean_orthologs"""

from divergence import resource_filename
from divergence.clean_orthologs import _extract_shared_orthologs, trim_and_concat_sicos, _align_sicos, \
_trim_alignments, _run_translatorx, _trim_alignment
from divergence.test import ECOLI_GENOMES, SALMONELLA_GENOMES, ECOLI_AND_SALMONELLA_GENOMES, fail_if_diff
import os
import unittest

class Test(unittest.TestCase):
    """Test class"""

    def test_cleanup(self):
        """For all 2 ecoli & 2 salmonella, create files per ortholog, assure cogs match, and run TranslatorX."""
        genomes = ECOLI_GENOMES[:2] + SALMONELLA_GENOMES[:2]
        groups_file = resource_filename(__name__, 'data/orthomcl/groups.txt')
        dna_files = [resource_filename(__name__, 'data/orthomcl/58191.ffn'),
                     resource_filename(__name__, 'data/orthomcl/58531.ffn'),
                     resource_filename(__name__, 'data/orthomcl/59245.ffn'),
                     resource_filename(__name__, 'data/orthomcl/59431.ffn')]

        trimmed_sicos, concatemers, stats_file = trim_and_concat_sicos(genomes, dna_files, groups_file)

        #Assert we get as much results as we expected
        self.failUnlessEqual(2586, len(trimmed_sicos))
        self.failUnlessEqual(len(genomes), len(concatemers))

        for trimmed_sico in trimmed_sicos:
            self.failUnless(0 < os.path.getsize(trimmed_sico))
        for concatemer in concatemers:
            self.failUnless(0 < os.path.getsize(concatemer))

        #Diff file results with expected results, rather than merely checking file existence
        expected_stats = resource_filename(__name__, 'data/orthomcl/stats.txt')
        fail_if_diff(self, expected_stats, stats_file)

    def test_trim_alignments(self):
        """Trim two alignments and compare them with the same alignments trimmed by hand."""
        sicos = [resource_filename(__name__, 'data/cleanup/sico_000000.ffn'),
                      resource_filename(__name__, 'data/cleanup/sico_000001.ffn')]

        #Align and trim
        alignments = _align_sicos(sicos)
        trimmed = _trim_alignments(alignments, '/dev/null')

        expected = [resource_filename(__name__, 'data/cleanup/sico_000000.nt_ali.fasta'),
                      resource_filename(__name__, 'data/cleanup/sico_000001.nt_ali.fasta')]

        #Compare using fail_if_diff, making use the exit code is 0 for no differences found
        fail_if_diff(self, expected[0], trimmed[0])
        fail_if_diff(self, expected[1], trimmed[1])

        #Run again, but now singular
        for sico, expected in zip(sicos, expected):
            alignment = _run_translatorx(sico)
            trimmed = _trim_alignment(alignment)[0]
            fail_if_diff(self, expected, trimmed)

    def test_filter_exact_set_10v10(self):
        """Test generating statistics for the exact set of genomes used for the groups.txt file."""
        groups_file = resource_filename(__name__, 'data/cleanup/groups-10-vs-10.txt')

        #Exact set of genomes that were used to generate the above groups.txt file
        genomes = ECOLI_AND_SALMONELLA_GENOMES

        print '10-vs-10 exact set'

        #Assertions
        sico, muco, non_shared = _extract_shared_orthologs(genomes, groups_file)
        self._check_values((2048, 64, 7389), sico, muco, non_shared)

    def test_filter_sub_set_10v10(self):
        """Test generating statistics for a sub set of genomes used for the groups.txt file."""
        groups_file = resource_filename(__name__, 'data/cleanup/groups-10-vs-10.txt')

        #Sub set of genomes that were used to generate the above groups.txt file
        genomes = ECOLI_GENOMES[:2] + SALMONELLA_GENOMES[:2]

        #Should probably give more hits than for 3v3 set, but still different  
        print '10-vs-10 sub set'

        #Assertions
        sico, muco, non_shared = _extract_shared_orthologs(genomes, groups_file)
        self._check_values((2586, 71, 7389), sico, muco, non_shared)

    def test_filter_exact_set_3v3(self):
        """Test generating statistics for the exact set of genomes used for the groups.txt file."""
        groups_file = resource_filename(__name__, 'data/cleanup/groups-3-vs-3.txt')

        #Exact set of genomes that were used to generate the above groups.txt file
        genomes = ECOLI_GENOMES[:3] + SALMONELLA_GENOMES[:3]

        print '3-vs-3 exact set'

        #Assertions
        sico, muco, non_shared = _extract_shared_orthologs(genomes, groups_file)
        self._check_values((2508, 77, 5340), sico, muco, non_shared)

    def test_filter_super_set_3v3(self):
        """Test generating statistics for a super set of genomes used for the groups.txt file."""
        groups_file = resource_filename(__name__, 'data/cleanup/groups-3-vs-3.txt')
        #Superset of genomes that were used to generate the above groups.txt file
        genomes = ECOLI_GENOMES[:4] + SALMONELLA_GENOMES[:4]

        print '3-vs-3 super set'

        #Assertions
        sico, muco, non_shared = _extract_shared_orthologs(genomes, groups_file)
        self._check_values((0, 0, 5340), sico, muco, non_shared)

    def test_filter_sub_set_3v3(self):
        """Test generating statistics for a sub set of genomes used for the groups.txt file."""
        groups_file = resource_filename(__name__, 'data/cleanup/groups-3-vs-3.txt')
        #Sub set of genomes that were used to generate the above groups.txt file
        genomes = ECOLI_GENOMES[:2] + SALMONELLA_GENOMES[:2]

        print '3-vs-3 sub set'

        #Assertions
        sico, muco, non_shared = _extract_shared_orthologs(genomes, groups_file)
        self._check_values((2593, 75, 5340), sico, muco, non_shared)

    def _check_values(self, (exp_sico, exp_muco, exp_total), sico, muco, non_shared):
        """Assert actual values match expected values. This method is here to reduce code duplication in each test."""
        total = len(sico) + len(muco) + len(non_shared)
        print len(sico), len(muco), len(non_shared), total
        self.failUnlessEqual(exp_sico, len(sico))
        self.failUnlessEqual(exp_muco, len(muco))
        self.failUnlessEqual(exp_total, len(sico) + len(muco) + len(non_shared))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_generate_statistics']
    unittest.main()
