#!/usr/local/bin/python2.7
# encoding: utf-8
'''
divergence.calculations_new -- calculate some values for aligned single copy orthologs

divergence.calculations_new is a module to calculate values for SICOs

@author:     Tim te Beek
@contact:    tim.te.beek@nbic.nl
'''

from __future__ import division
from Bio.Data import CodonTable
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import Counter, defaultdict
from divergence import CODON_TABLE_ID
from itertools import product
import os
import sys



__all__ = []
__version__ = 0.1
__date__ = '2013-01-16'
__updated__ = '2013-01-16'

DEBUG = 1

# Premise
# - Some duplication is OK if it helps clarity
# - Do not repeatedly pass around the same arguments
# - Use constants for calculated value keys
# - Separate each calculation run from all others
# - Argument types should always be clear


# Field name constants for clearer usage patterns
ORTHOLOG = ''
PRODUCT = ''
COG_DIGITS = ''
COG_LETTERS = ''
CODONS = ''
NON_SYNONYMOUS_SITES = ''
NON_SYNONYMOUS_POLYMORPHISMS = ''
NON_SYNONYMOUS_SFS_SINGLETONS = ''
NON_SYNONYMOUS_SFS_DOUBLETONS = ''
NON_SYNONYMOUS_SFS_TRIPLETONS = ''
SYNONYMOUS_SITES = ''
SYNONYMOUS_POLYMORPHISMS = ''
SYNONYMOUS_SFS_SINGLETONS = ''
SYNONYMOUS_SFS_DOUBLETONS = ''
SYNONYMOUS_SFS_TRIPLETONS = ''
FOUR_FOLD_SYNONYMOUS_SITES = ''
FOUR_FOLD_SYNONYMOUS_POLYMORPHISMS = ''
FOUR_FOLD_SYNONYMOUS_SFS_SINGLETONS = ''
FOUR_FOLD_SYNONYMOUS_SFS_DOUBLETONS = ''
FOUR_FOLD_SYNONYMOUS_SFS_TRIPLETONS = ''
MULTIPLE_SITE_POLYMORPHISMS = ''
COMPLEX_CODONS = ''
DN = ''
DS = ''
PHIPACK_SITES = ''
PHI = ''
MAX_CHI_2 = ''
NSS = ''
PI = ''
PI_NONSYN = ''
PI_SYN = ''
PI_4_FOLD_SYN = ''
THETA = ''
NEUTRALITY_INDEX = ''
DOS = ''


def _calc_pi(nr_of_strains, sequence_lengths, site_freq_spec):
    """
    New, improved: n/(n-1) * Sum( Pj * 2 * j/n * (1-j/n), {i,1,Floor((n-1)/2)})

    Pi: n/(n-1) * Sum[ D(i) * 2 i/n (1-i/n), i from 1 to RoundDown(n/2)]
    where n is number of strains
    and D(i) is the number of polymorphisms present in i of n strains
    finally divide everything by the number of sites
    """
    print site_freq_spec
    return (nr_of_strains
                     / (nr_of_strains - 1)
                     * sum(site_freq_spec.get(i, 0)
                           * 2 * i / nr_of_strains
                           * (1 - i / nr_of_strains)
                           for i in range(1, (nr_of_strains - 1) // 2 + 1))  # +1 as range excludes stop value
                     ) / sequence_lengths

def _site_freq_spec(clade_calcs):
    '''Add the full site frequency spectrum to a clade_calcs instance.'''
    site_freq_spec = defaultdict(int)

    # Loop over every position in alignment to determine the bases present
    for pos in range(clade_calcs.sequence_lengths):
        counts = Counter(clade_calcs.alignment[:, pos])

        # If we only find a single base, we can ignore this site
        if len(counts) == 1:
            continue

        # Remove the most prevalent base from the counts
        most_prevalent = max(counts.keys(), key=lambda x: counts[x])
        del counts[most_prevalent]

        # Update SFS with the remaining bases
        for occurrences in counts.values():
            site_freq_spec[occurrences] += 1

    clade_calcs.values[PI] = _calc_pi(clade_calcs.nr_of_strains, clade_calcs.sequence_lengths, site_freq_spec)

# Using the standard NCBI Bacterial, Archaeal and Plant Plastid Code translation table (11)
BACTERIAL_CODON_TABLE = CodonTable.unambiguous_dna_by_id.get(CODON_TABLE_ID)

def _four_fold_degenerate_patterns():
    """Find patterns of 4-fold degenerate codons, wherein all third site substitutions code for the same amino acid."""
    letters = BACTERIAL_CODON_TABLE.nucleotide_alphabet.letters
    # Any combination of letters of length two
    for site12 in [''.join(prod) for prod in product(letters, repeat=2)]:
        # 4-fold when the length of the unique encoded amino acids for all possible third site nucleotides is exactly 1
        if 1 == len(set([BACTERIAL_CODON_TABLE.forward_table.get(site12 + site3) for site3 in letters])):
            # Add regular expression pattern to the set of patterns
            yield '{0}[{1}]'.format(site12, letters)

FOUR_FOLD_DEGENERATE_PATTERN = '|'.join(_four_fold_degenerate_patterns())

def _codon_site_freq_spec(clade_calcs):
    ''''''
    multiple_site_polymorphisms = 0
    mixed_synonymous_polymorphisms = 0
    
    global_sfs = defaultdict(int)

    synonymous_sfs = defaultdict(int)
    non_synonymous_sfs = defaultdict(int)
    four_fold_syn_sfs = defaultdict(int)

    # Calculate sequence_lengths here so we can handle alignments that are not multiples of three
    sequence_lengths = clade_calcs.sequence_lengths - clade_calcs.sequence_lengths % 3

    # Split into codon_alignments
    codon_alignments = (clade_calcs.alignment[:, index:index + 3] for index in range(0, sequence_lengths, 3))
    for codon_alignment in codon_alignments:
        # Get string representations of codons for simplicity
        codons = [str(seqr.seq) for seqr in codon_alignment]

        # Skip when all codons are the same
        if len(set(codons)) == 1:
            continue

        # As per AEW: Skip codons with gaps, and codons with unresolved bases: Basically anything but ACGT
        if 0 < len(''.join(codons).translate(None, 'ACGTactg')):
            continue

        # Skip codons where any of the alignment codons is a stopcodon, same as in codeml
        for codon in codons:
            if codon in BACTERIAL_CODON_TABLE.stop_codons:
                continue

        # Determine variation per site
        per_site_usage = [Counter(codon_alignment[:, site]) for site in range(3)]

        # Determine which sites contain polymorphisms
        polymorph_site_usages = [usage for usage in per_site_usage if 1 < len(usage)]

        # Skip codons where multiple sites contain polymorphisms
        if 1 < len(polymorph_site_usages):
            multiple_site_polymorphisms += 1
            continue

        # Extract the polymorphic site
        polymorph_site_usage = polymorph_site_usages[0]

        # Find the most prevalent base from the counts so we can ignore it for the SFS
        most_prevalent_base = max(polymorph_site_usage.keys(), key=lambda x: polymorph_site_usage[x])

        # Calculate the local site frequency spectrum
        local_sfs = defaultdict(int)
        for base, counts in polymorph_site_usage.items():
            if base == most_prevalent_base:
                continue
            else:
                local_sfs[counts] += 1

        def add_dict_to_dict(target, source):
            '''Add values from source to target'''
            for key, value in source.iteritems():
                target[key] += value

        # Global SFS takes it values from the local SFS, no further filtering applied
        add_dict_to_dict(global_sfs, local_sfs)

        # Retrieve translations of codons now that inconclusive & stop-codons have been removed
        translations = Counter(BACTERIAL_CODON_TABLE.forward_table.get(codon) for codon in codons)

        if len(translations) == 1:
            # All mutations are synonymous
            add_dict_to_dict(synonymous_sfs, local_sfs)
        else:
            if len(translations) == len(polymorph_site_usage):
                # Multiple translations, one per change in base
                add_dict_to_dict(non_synonymous_sfs, local_sfs)
            else:
                # Number of translations & number of different bases do not match: Both syn and non syn changes found
                mixed_synonymous_polymorphisms += 1


    print 'Global'
    print _calc_pi(clade_calcs.nr_of_strains, clade_calcs.sequence_lengths, global_sfs)

    print '\nSyn'
    print _calc_pi(clade_calcs.nr_of_strains, clade_calcs.sequence_lengths, synonymous_sfs)

    print '\nNon Syn'
    print _calc_pi(clade_calcs.nr_of_strains, clade_calcs.sequence_lengths, non_synonymous_sfs)

    print '\nMultiple site poly', multiple_site_polymorphisms


class clade_calcs(object):
    '''Perform the calculations specific a single clade.'''

    alignment = None
    nr_of_strains = None
    sequence_lengths = None

    values = {}

    def __init__(self, alignment):
        self.alignment = alignment
        self.nr_of_strains = len(alignment)
        self.sequence_lengths = len(alignment[0])


def main(argv=None):  # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_usage = '''%s

  Created by Tim te Beek on %s.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_usage, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %(default)s]")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        parser.add_argument(dest="paths", help="paths to folder(s) with source file(s) [default: %(default)s]", metavar="path", nargs='+')

        # Process arguments
        args = parser.parse_args()

        paths = args.paths
        verbose = args.verbose

        if verbose > 0:
            print("Verbose mode on")

        for inpath in paths:
            ### do something with inpath ###
            print(inpath)
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        if DEBUG:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
        sys.argv.append("-r")
    sys.exit(main())
