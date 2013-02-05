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
from divergence import CODON_TABLE_ID, find_cogs_in_sequence_records
from itertools import product
import logging
import os
import re
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
ORTHOLOG = 'ortholog'
PRODUCT = 'product'
COG_DIGITS = 'cog digits'
COG_LETTERS = 'cog letters'
CODONS = 'codons'
NON_SYNONYMOUS_SITES = 'non-synonymous sites'
NON_SYNONYMOUS_POLYMORPHISMS = 'non-synonymous polymorphisms'
SYNONYMOUS_SITES = 'synonymous sites'
SYNONYMOUS_POLYMORPHISMS = 'synonymous polymorphisms'
FOUR_FOLD_SYNONYMOUS_SITES = '4-fold synonymous sites'
FOUR_FOLD_SYNONYMOUS_POLYMORPHISMS = '4-fold synonymous polymorphisms'
MULTIPLE_SITE_POLYMORPHISMS = 'multiple site polymorphisms'
COMPLEX_CODONS = 'complex codons (with both synonymous and non-synonymous polymorphisms segregating)'
DN = 'Dn'
DS = 'Ds'
PHIPACK_SITES = 'PhiPack sites'
PHI = 'Phi'
MAX_CHI_2 = 'Max Chi^2'
NSS = 'NSS'
PI = 'Pi'
PI_NONSYN = 'Pi nonsyn'
PI_SYN = 'Pi syn'
PI_4_FOLD_SYN = 'Pi 4-fold syn'
THETA = 'Theta'
NEUTRALITY_INDEX = 'neutrality index'
DOS = 'DoS'


def _extract_cog_digits_and_letters(clade_calcs):
    '''Add the COG digits and letters to the clade_calcs.values dictionary for all strains in clade_calcs.alignment.'''
    cog_digits = []
    cog_letters = []
    for cog in find_cogs_in_sequence_records(clade_calcs.alignment):
        # Match digits and letters separately
        matchobj = re.match('(COG[0-9]+)([A-Z]*)', cog)
        if matchobj:
            cog_digits.append(matchobj.groups()[0])
            cog_letters.append(matchobj.groups()[1])
    # Join the found digits and letters using a comma
    clade_calcs.values[COG_DIGITS] = ','.join(cog_digits)
    clade_calcs.values[COG_LETTERS] = ','.join(cog_letters)


def _calc_pi(nr_of_strains, sequence_lengths, site_freq_spec):
    """
    New, improved: n/(n-1) * Sum( Pj * 2 * j/n * (1-j/n), {i,1,Floor((n-1)/2)})

    Pi: n/(n-1) * Sum[ D(i) * 2 i/n (1-i/n), i from 1 to RoundDown(n/2)]
    where n is number of strains
    and D(i) is the number of polymorphisms present in i of n strains
    finally divide everything by the number of sites
    """
    pi = (nr_of_strains
                     / (nr_of_strains - 1)
                     * sum(site_freq_spec.get(i, 0)
                           * 2 * i / nr_of_strains
                           * (1 - i / nr_of_strains)
                           for i in range(1, (nr_of_strains - 1) // 2 + 1))  # +1 as range excludes stop value
                     ) / sequence_lengths
    logging.debug('Arguments:\n\tnr of strains: %s,\n\tsequence lengths: %s,\n\tsfs: %s',
                  nr_of_strains,
                  sequence_lengths,
                  site_freq_spec)
    logging.debug('Gave a Pi value of: %s', pi)
    return pi

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
    '''Find patterns of 4-fold degenerate codons, wherein all third site substitutions code for the same amino acid.'''
    letters = BACTERIAL_CODON_TABLE.nucleotide_alphabet.letters
    # Any combination of letters of length two
    for site12 in [''.join(prod) for prod in product(letters, repeat=2)]:
        # 4-fold when the length of the unique encoded amino acids for all possible third site nucleotides is exactly 1
        if 1 == len(set([BACTERIAL_CODON_TABLE.forward_table.get(site12 + site3) for site3 in letters])):
            # Add regular expression pattern to the set of patterns
            yield site12

FOUR_FOLD_DEGENERATE_PATTERN = '({onetwo})[{three}]'.format(onetwo='|'.join(sorted(_four_fold_degenerate_patterns())),
                                                            three=BACTERIAL_CODON_TABLE.nucleotide_alphabet.letters)

def _codon_site_freq_spec(clade_calcs):
    '''Site frequency spectrum calculations for full, syn, non-syn and 4-fold syn sites.'''
    four_fold_synonymous_sites = 0
    multiple_site_polymorphisms = 0
    mixed_synonymous_polymorphisms = 0

    stop_codons = 0
    codons_with_unresolved_bases = 0

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
            # Increase number of four fold synonymous sites if the codons match; No SFS to add as all codons are equal
            if re.match(FOUR_FOLD_DEGENERATE_PATTERN, codons[0]):
                four_fold_synonymous_sites += 1
            continue

        # As per AEW: Skip codons with gaps, and codons with unresolved bases: Basically anything but ACGT
        if 0 < len(''.join(codons).translate(None, 'ACGTactg')):
            codons_with_unresolved_bases += 1
            continue

        # Skip codons where any of the alignment codons is a stopcodon, same as in codeml
        for codon in codons:
            if codon in BACTERIAL_CODON_TABLE.stop_codons:
                stop_codons += 1
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

            # Check if these codons also match the four fold synonymous pattern
            if all(re.match(FOUR_FOLD_DEGENERATE_PATTERN, codon) for codon in codons):
                add_dict_to_dict(four_fold_syn_sfs, local_sfs)
                four_fold_synonymous_sites += 1
        else:
            if len(translations) == len(polymorph_site_usage):
                # Multiple translations, one per change in base
                add_dict_to_dict(non_synonymous_sfs, local_sfs)
            else:
                # Number of translations & number of different bases do not match: Both syn and non syn changes found
                mixed_synonymous_polymorphisms += 1

        # Implicitly continue with next iteration

    # Add calculations to values dictionary
    logging.debug('%s SFS: %s', PI, global_sfs)
    clade_calcs.values[PI] = _calc_pi(clade_calcs.nr_of_strains, clade_calcs.sequence_lengths, global_sfs)

    logging.debug('%s SFS: %s', PI_SYN, synonymous_sfs)
    clade_calcs.values[PI_SYN] = _calc_pi(clade_calcs.nr_of_strains, clade_calcs.sequence_lengths, synonymous_sfs)
    clade_calcs.values[SYNONYMOUS_POLYMORPHISMS] = sum(synonymous_sfs.values())

    logging.debug('%s SFS: %s', PI_NONSYN, non_synonymous_sfs)
    clade_calcs.values[PI_NONSYN] = _calc_pi(clade_calcs.nr_of_strains, clade_calcs.sequence_lengths, non_synonymous_sfs)
    clade_calcs.values[NON_SYNONYMOUS_POLYMORPHISMS] = sum(non_synonymous_sfs.values())

    logging.debug('%s SFS: %s', PI_4_FOLD_SYN, four_fold_syn_sfs)
    clade_calcs.values[PI_4_FOLD_SYN] = _calc_pi(clade_calcs.nr_of_strains, clade_calcs.sequence_lengths, four_fold_syn_sfs)
    clade_calcs.values[FOUR_FOLD_SYNONYMOUS_POLYMORPHISMS] = sum(four_fold_syn_sfs.values())

    # TODO Also add the SINGLETON, DOUBLETON, TRIPLETON, etc values here

    # Add tallies to values dictionary
    clade_calcs.values[FOUR_FOLD_SYNONYMOUS_SITES] = four_fold_synonymous_sites
    clade_calcs.values[MULTIPLE_SITE_POLYMORPHISMS] = multiple_site_polymorphisms
    clade_calcs.values[COMPLEX_CODONS] = mixed_synonymous_polymorphisms

    # Log debug statistics
    logging.debug('stop_codons: %s', stop_codons)
    logging.debug('codons_with_unresolved_bases: %s', codons_with_unresolved_bases)


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

        # The most basic calculation added to the output file
        clade_calcs.values[CODONS] = self.sequence_lengths // 3

        # TODO extract product for most recent gene


def run_calculations(genomes_a_file,
                     genomes_b_file,
                     sicozip_file,
                     table_a_dest,
                     table_b_dest,
                     append_odd_even=False):
    ''''''
    # TODO All stubs below
    # parse genomes in genomes_x_files

    # extract ortholog files from sicozip
    # loop over orthologs
    # # calculate phipack values for combined aligments
    # # split alignments
    # # # calculate codeml values
    # # # COG digits and letters
    # # # SFS related values
    # # # additional calculations: theta, ni, DoS
    # # # bootstrapping NI

    # odd even tables




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

        # Arguments specific to calculations
        parser.add_argument('genomes_a')
        parser.add_argument('genomes_b')
        parser.add_argument('sicozip')
        parser.add_argument('table_a')
        parser.add_argument('table_b')

        parser.add_argument('-a', '--append-odd-even', action='store_true',
                            help='append separate tables calculated for odd and even codons of ortholog alignments (default: False)')

        # Process arguments
        args = parser.parse_args()

        if args.verbose > 0:
            print("Verbose mode on")
            logging.root.setLevel(logging.DEBUG)

        run_calculations(args.genomes_a,
                         args.genomes_b,
                         args.sicozip,
                         args.table_a,
                         args.table_b,
                         args.append_odd_even)

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
    sys.exit(main())
