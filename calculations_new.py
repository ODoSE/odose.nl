#!/usr/local/bin/python2.7
# encoding: utf-8
'''
divergence.calculations_new -- calculate some values for aligned single copy orthologs

divergence.calculations_new is a module to calculate values for SICOs

@author:     Tim te Beek
@contact:    tim.te.beek@nbic.nl
'''

from __future__ import division
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable
from argparse import ArgumentParser, RawDescriptionHelpFormatter, ArgumentTypeError
from collections import Counter, defaultdict
from divergence import CODON_TABLE_ID, find_cogs_in_sequence_records, get_most_recent_gene_name, \
    extract_archive_of_files, create_directory
from divergence.run_codeml import run_codeml, parse_codeml_output
from divergence.run_phipack import run_phipack
from divergence.select_taxa import select_genomes_by_ids
from itertools import product
from numpy import mean
from random import choice
import logging
import os
import re
import shutil
import sys
import tempfile


__all__ = []
__version__ = 0.1
__date__ = '2013-01-16'
__updated__ = '2013-02-11'

DEBUG = 0

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
PI = 'Pi'
GLOBAL_SFS = 'global sfs'
NON_SYNONYMOUS_PI = 'Pi nonsyn'
NON_SYNONYMOUS_SFS = 'non-synonymous sfs'
NON_SYNONYMOUS_SITES = 'non-synonymous sites'
NON_SYNONYMOUS_POLYMORPHISMS = 'non-synonymous polymorphisms'
SYNONYMOUS_PI = 'Pi syn'
SYNONYMOUS_SFS = 'synonymous sfs'
SYNONYMOUS_SITES = 'synonymous sites'
SYNONYMOUS_POLYMORPHISMS = 'synonymous polymorphisms'
FOUR_FOLD_SYNONYMOUS_PI = 'Pi 4-fold syn'
FOUR_FOLD_SYNONYMOUS_SFS = '4-fold synonymous sfs'
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
THETA = 'Theta'
DS_PN_PS_DS = 'Ds*Pn/(Ps+Ds)'
DN_PS_PS_DS = 'Dn*Ps/(Ps+Ds)'
NEUTRALITY_INDEX = 'neutrality index'
DOS = 'DoS'


def _get_nton_name(nton, prefix=''):
    """Given the number of strains in which a polymorphism/substitution is found, give the appropriate SFS name."""
    named = {1: 'single', 2: 'double', 3: 'triple', 4: 'quadruple', 5: 'quintuple'}
    middle = named.get(nton, str(nton) + '-')
    return prefix + middle + 'tons'


def _get_column_headers(max_nton):
    '''Get the column headers in the order they need to appear in the output data file.'''

    def _write_sfs_column_names(prefix):
        '''Return named columns for SFS singleton, doubleton, tripleton, etc.. upto max_nton.'''
        for number in range(1, max_nton + 1):
            yield _get_nton_name(number, prefix)

    headers = [
               # Ortholog is first
               ORTHOLOG,

               # Gene name is captured as product of the orthologous DNA sequence
               PRODUCT,

               # Split COG into first part with numbers and second part with letters
               COG_DIGITS, COG_LETTERS,
               # Some initial values
               CODONS]

    # Non_synonymous
    headers.extend([NON_SYNONYMOUS_SITES, NON_SYNONYMOUS_POLYMORPHISMS])
    headers.extend(_write_sfs_column_names(NON_SYNONYMOUS_SFS + ' '))

    # synonymous
    headers.extend([SYNONYMOUS_SITES, SYNONYMOUS_POLYMORPHISMS])
    headers.extend(_write_sfs_column_names(SYNONYMOUS_SFS + ' '))

    # 4-fold synonymous
    headers.extend([FOUR_FOLD_SYNONYMOUS_SITES, FOUR_FOLD_SYNONYMOUS_POLYMORPHISMS])
    headers.extend(_write_sfs_column_names(FOUR_FOLD_SYNONYMOUS_SFS + ' '))

    headers.extend([
                    # Miscellaneous additional statistics
                    MULTIPLE_SITE_POLYMORPHISMS,
                    COMPLEX_CODONS,
                    # PAML
                    DN, DS,
                    # PhiPack
                    PHIPACK_SITES, PHI, MAX_CHI_2, NSS,

                    # Measures of nucleotide diversity per SICO
                    PI,
                    NON_SYNONYMOUS_PI,
                    SYNONYMOUS_PI,
                    FOUR_FOLD_SYNONYMOUS_PI,
                    THETA,

                    # Final _two_ values here are ignored when calculation sums over complete table
                    NEUTRALITY_INDEX,
                    DOS])
    return headers


def _write_intro_to_file(table_a_dest):
    '''Explain the presence of three tables in one file if we're also calculating odd & even codons'''
    with open(table_a_dest, 'w') as write_handle:
        write_handle.write('''#First table contains calculations for all codon
#Second table contains calculations for odd codons only
#Third table contains calculations for even codons only''')


def _write_to_file(table_a_dest,
                   genome_ids_a,
                   genome_ids_b,
                   common_prefix_a,
                   common_prefix_b,
                   calculations):
    '''

    :param table_a_dest:
    :type table_a_dest: filename
    :param genome_ids_a:
    :type genome_ids_a: list or genome ids
    :param common_prefix_a:
    :type common_prefix_a: string
    :param common_prefix_b:
    :type common_prefix_b: string
    :param calculations:
    :type calculations: list of clade_calcs instances
    '''
    with open(table_a_dest, 'a') as write_handle:
        # Print introduction about the strain comparison
        write_handle.write('#{} {} strains compared with {} {} strains\n'.format(len(genome_ids_a),
                                                                                 common_prefix_a,
                                                                                 len(genome_ids_b),
                                                                                 common_prefix_b))
        # Print the genome IDs involved in each of the strains
        write_handle.write('#IDs {}: {}\n'.format(common_prefix_a,
                                                  ', '.join(genome_ids_a)))
        write_handle.write('#IDs {}: {}\n'.format(common_prefix_b,
                                                  ', '.join(genome_ids_b)))

        # Print column headers for the data to come
        max_nton = len(genome_ids_a) // 2
        headers = _get_column_headers(max_nton)
        write_handle.write('#' + '\t'.join(headers))
        write_handle.write('\n')

        # Print data rows
        format_str = '\t'.join('{{{}}}'.format(key) for key in headers)
        from string import Formatter
        formatter = Formatter()
        for clade_calcs in calculations:
            write_handle.write(formatter.vformat(format_str, None, clade_calcs.values))
            write_handle.write('\n')


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


def _get_codeml_values(alignment_a, alignment_b):
    '''Get the codeml values for running the first sequences of both alignment a & b through codeml and return dict.'''
    # Run codeml to calculate values for dn & ds
    subdir = tempfile.mkdtemp(prefix='codeml_')
    codeml_file = run_codeml(subdir, alignment_a, alignment_b)
    codeml_values_dict = parse_codeml_output(codeml_file)
    shutil.rmtree(subdir)

    # convert poorly legible keys to better ones
    codeml_values_dict[NON_SYNONYMOUS_SITES] = codeml_values_dict['N']
    codeml_values_dict[SYNONYMOUS_SITES] = codeml_values_dict['S']

    return codeml_values_dict


def _calc_pi(nr_of_strains, nr_of_sites, site_freq_spec):
    """
    New, improved: n/(n-1) * Sum( Pj * 2 * j/n * (1-j/n), {i,1,Floor((n-1)/2)})

    Pi: n/(n-1) * Sum[ D(i) * 2 i/n (1-i/n), i from 1 to RoundDown(n/2)]
    where n is number of strains
    and D(i) is the number of polymorphisms present in i of n strains
    finally divide everything by the number of sites
    """
    if not site_freq_spec:
        return 0.0
    logging.debug('args:\tstrains: %s,\tseqlen: %s,\tsfs: %s',
                  nr_of_strains,
                  nr_of_sites,
                  site_freq_spec)
    pi = (nr_of_strains
                     / (nr_of_strains - 1)
                     * sum(site_freq_spec.get(i, 0)
                           * 2 * i / nr_of_strains
                           * (1 - i / nr_of_strains)
                           for i in range(1, (nr_of_strains - 1) // 2 + 1))  # +1 as range excludes stop value
                     ) / nr_of_sites
    logging.debug('Gave a Pi value of: %s', pi)
    return pi


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
                four_fold_synonymous_sites += 1
                add_dict_to_dict(four_fold_syn_sfs, local_sfs)
        else:
            if len(translations) == len(polymorph_site_usage):
                # Multiple translations, one per change in base
                add_dict_to_dict(non_synonymous_sfs, local_sfs)
            else:
                # Number of translations & number of different bases do not match: Both syn and non syn changes found
                mixed_synonymous_polymorphisms += 1

        # Implicitly continue with next iteration

    # Add SFS & Pi calculations to values dictionary
    # Synonymous
    clade_calcs.values[GLOBAL_SFS] = global_sfs
    clade_calcs.values[PI] = _calc_pi(clade_calcs.nr_of_strains, clade_calcs.sequence_lengths, global_sfs)

    clade_calcs.values[SYNONYMOUS_SFS] = synonymous_sfs
    clade_calcs.values[SYNONYMOUS_POLYMORPHISMS] = sum(synonymous_sfs.values())
    clade_calcs.values[SYNONYMOUS_PI] = _calc_pi(clade_calcs.nr_of_strains,
                                                 clade_calcs.values[SYNONYMOUS_SITES],
                                                 synonymous_sfs)
    for nton, value in synonymous_sfs.items():
        clade_calcs.values[_get_nton_name(nton, SYNONYMOUS_SFS + ' ')] = value

    # Non synonymous
    clade_calcs.values[NON_SYNONYMOUS_SFS] = non_synonymous_sfs
    clade_calcs.values[NON_SYNONYMOUS_POLYMORPHISMS] = sum(non_synonymous_sfs.values())
    clade_calcs.values[NON_SYNONYMOUS_PI] = _calc_pi(clade_calcs.nr_of_strains,
                                                     clade_calcs.values[NON_SYNONYMOUS_SITES],
                                                     non_synonymous_sfs)
    for nton, value in non_synonymous_sfs.items():
        clade_calcs.values[_get_nton_name(nton, NON_SYNONYMOUS_SFS + ' ')] = value

    # 4-fold synonymous
    clade_calcs.values[FOUR_FOLD_SYNONYMOUS_SFS] = four_fold_syn_sfs
    clade_calcs.values[FOUR_FOLD_SYNONYMOUS_POLYMORPHISMS] = sum(four_fold_syn_sfs.values())
    clade_calcs.values[FOUR_FOLD_SYNONYMOUS_PI] = _calc_pi(clade_calcs.nr_of_strains,
                                                           four_fold_synonymous_sites,
                                                           four_fold_syn_sfs)
    for nton, value in four_fold_syn_sfs.items():
        clade_calcs.values[_get_nton_name(nton, FOUR_FOLD_SYNONYMOUS_SFS + ' ')] = value


    # Add tallies to values dictionary
    clade_calcs.values[FOUR_FOLD_SYNONYMOUS_SITES] = four_fold_synonymous_sites
    clade_calcs.values[MULTIPLE_SITE_POLYMORPHISMS] = multiple_site_polymorphisms
    clade_calcs.values[COMPLEX_CODONS] = mixed_synonymous_polymorphisms

    # Log debug statistics
    if stop_codons:
        logging.debug('stop_codons: %s', stop_codons)
    if codons_with_unresolved_bases:
        logging.debug('codons_with_unresolved_bases: %s', codons_with_unresolved_bases)


def _extract_genome_ids_and_common_prefix(genomes_file):
    '''From a genome ids file extract the genome ids and the common name prefix for all genomes.'''
    with open(genomes_file) as read_handle:
        lines = [line.strip() for line in read_handle]
        genome_ids = [line.split()[0] for line in lines]
        common_prefix = os.path.commonprefix([line.split('\t')[1] for line in lines]).strip()
        return genome_ids, common_prefix


def _add_combined_calculations(clade_calcs):
    '''Add additional deduced calculations: theta, ni, DoS...'''

    # 16. Direction of selection = Dn/(Dn+Ds) - Pn/(Pn+Ps)
    paml_total_substitutions = clade_calcs.values[DN] + clade_calcs.values[DS]
    total_polymorphisms = clade_calcs.values[NON_SYNONYMOUS_POLYMORPHISMS] + clade_calcs.values[SYNONYMOUS_POLYMORPHISMS]
    # Prevent divide by zero by checking both values above are not null
    if paml_total_substitutions != 0 and total_polymorphisms != 0:
        clade_calcs.values[DOS] = (clade_calcs.values[DN] / paml_total_substitutions
                                                 - clade_calcs.values[NON_SYNONYMOUS_POLYMORPHISMS] / total_polymorphisms)
    else:
        # "the direction of selection is undefined if either Dn+Ds or Pn+Ps are zero": None or Not a Number?
        clade_calcs.values[DOS] = None

    # Theta
    # Watterson's estimator of theta: S / (L * harmonic)
    # where the harmonic is Sum[ 1 / i, i from 1 to n - 1 ]
    harmonic = sum(1 / i for i in range(1, clade_calcs.nr_of_strains))
    clade_calcs.values[THETA] = sum(clade_calcs.values[GLOBAL_SFS].values()) / (clade_calcs.sequence_lengths
                                                                                * harmonic)

    # NI (parts)
    # These values will end up contributing to the Neutrality Index through NI = Sum(X) / Sum(Y)
    ps_plus_ds = (clade_calcs.values[SYNONYMOUS_POLYMORPHISMS] + clade_calcs.values[DS])
    if ps_plus_ds:
        # X = Ds*Pn/(Ps+Ds)
        clade_calcs.values[DS_PN_PS_DS] = (clade_calcs.values[DS]
                                               * clade_calcs.values[NON_SYNONYMOUS_POLYMORPHISMS]
                                               / ps_plus_ds)
        # Y = Dn*Ps/(Ps+Ds)
        clade_calcs.values[DN_PS_PS_DS] = (clade_calcs.values[DN]
                                               * clade_calcs.values[SYNONYMOUS_POLYMORPHISMS]
                                               / ps_plus_ds)
    else:
        clade_calcs.values[DS_PN_PS_DS] = None
        clade_calcs.values[DN_PS_PS_DS] = None


class Statistic(object):
    '''Helper class for general statistics that mimics the characterics of clade_calcs.'''
    def __init__(self, name):
        '''The name argument is stored as ORTHOLOG in the values defaultdict, so it shows up correctly in output.'''
        self.values = defaultdict(lambda:'')
        self.values[ORTHOLOG] = name


def _calculcate_mean_and_averages(calculations, max_nton):
    '''Calculate the sum and mean data rows for a subset of numerical data columns, and append them to calculations.'''
    headers = _get_column_headers(max_nton)

    # calculate the sum for a subset of headers
    sum_stats = Statistic('sum')
    for header in headers[5:-2]:
        sum_stats.values[header] = sum(clade_calcs.values[header]
                                       for clade_calcs in calculations
                                       if clade_calcs.values[header] != None)

    # calculate the average for a subset of headers
    mean_stats = Statistic('mean')
    for header in headers[5:-2] + [DOS]:
        # XXX fromnumeric.py:37: RuntimeWarning: invalid value encountered in double_scalars
        mean_stats.values[header] = mean([clade_calcs.values[header]
                                          for clade_calcs in calculations
                                          if clade_calcs.values[header] != None])

    # append statistics now that we're no longer looping over them
    return sum_stats, mean_stats


def _bootstrap(sum_dspn, sum_dnps):
    """Bootstrap by gene to get to confidence scores for Neutrality Index."""
    def _sample_with_replacement(sample_set, sample_size=None):
        """Sample sample_size items from sample_set, or len(sample_set) items if sample_size is None (default)."""
        if sample_size is None:
            sample_size = len(sample_set)
        samples = []
        while len(samples) < sample_size:
            samples.append(choice(sample_set))
        return samples

    # "to get the confident interval on this you need to boostrap by gene - i.e. if we have 1000 genes, we form a
    # boostrap sample by resampling, with replacement 1000 genes from the original sample; recalculate NI and repeat
    # 1000 times; the SE on the estimate is the standard deviation across bootstraps, and your 95% confidence
    # interval van be obtained by sorting the values and taking the 25t and 975th values"
    ni_values = []
    while len(ni_values) < len(sum_dspn):
        ni_values.append(sum(_sample_with_replacement(sum_dspn)) /
                         sum(_sample_with_replacement(sum_dnps)))

    # 95 percent of values fall between n*.025th element & n*.975th element when NI values are sorted
    ni_values = sorted(ni_values)
    lower_limit = int(round(0.025 * (len(ni_values) - 1)))
    upper_limit = int(round(0.975 * (len(ni_values) - 1)))
    return ni_values[lower_limit], ni_values[upper_limit]


def _neutrality_indices(calculations):
    '''Return the statistics for Neutrality index. It adds the actual value, and two bootstrapped 95% values.'''
    # Neutrality Index = Sum(X = Ds*Pn/(Ps+Ds)) / Sum(Y = Dn*Ps/(Ps+Ds))

    x_values = [clade_calcs.values[DS_PN_PS_DS]
                for clade_calcs in calculations
                if clade_calcs.values[DS_PN_PS_DS] != None]
    y_values = [clade_calcs.values[DN_PS_PS_DS]
                for clade_calcs in calculations
                if clade_calcs.values[DN_PS_PS_DS] != None]

    sum_x = sum(x_values)
    sum_y = sum(y_values)

    ni_stats = Statistic('NI')
    ni_lower_stats = Statistic('NI 95% lower limit')
    ni_upper_stats = Statistic('NI 95% upper limit')

    if sum_y:
        ni_stats.values[NEUTRALITY_INDEX] = sum_x / sum_y

        # Find lower and upper limits within which 95% of values fall, by using bootstrapping statistics
        lower_95perc_limit, upper_95perc_limit = _bootstrap(x_values, y_values)
        ni_lower_stats.values[NEUTRALITY_INDEX] = lower_95perc_limit
        ni_upper_stats.values[NEUTRALITY_INDEX] = upper_95perc_limit
    else:
        # We could not calculate Neutrality index because Sum(Y = Dn*Ps/(Ps+Ds)) was zero
        msg = 'Failed to calculate Neutrality index because divisor Sum(Dn*Ps/(Ps+Ds)) was zero'
        logging.warn(msg)
        ni_stats.values[NEUTRALITY_INDEX] = msg

    return ni_stats, ni_lower_stats, ni_upper_stats


class clade_calcs(object):
    '''Perform the calculations specific a single clade.'''

    alignment = None
    nr_of_strains = None
    sequence_lengths = None

    values = None

    def __init__(self, alignment, genomes):
        self.alignment = alignment
        self.nr_of_strains = len(alignment)
        self.sequence_lengths = len(alignment[0])

        self.values = defaultdict(int)

        # The most basic calculation added to the output file
        self.values[CODONS] = self.sequence_lengths // 3

        # Get the most recent gene name for the strains in a given clade_calcs instance
        self.values[PRODUCT] = get_most_recent_gene_name(genomes, self.alignment)


def _table_calculations(genome_ids_a, genome_ids_b, sico_files, phipack_values):
    '''Perform calculations for comparsion of genome_ids_a with genome_ids_b.'''
    # retrieve genomes once for both
    genomes_a = select_genomes_by_ids(genome_ids_a).values()

    # dictionary to hold the values calculated per file
    calculations = []
    # loop over orthologs
    for sico_file in sico_files:
        # parse alignment
        alignment = AlignIO.read(sico_file, 'fasta')

        # split alignments
        alignment_a = MultipleSeqAlignment(seqr for seqr in alignment if seqr.id.split('|')[0] in genome_ids_a)
        alignment_b = MultipleSeqAlignment(seqr for seqr in alignment if seqr.id.split('|')[0] in genome_ids_b)

        # calculate codeml values
        codeml_values = _get_codeml_values(alignment_a, alignment_b)

        # create gathering instance of clade_calcs
        instance = clade_calcs(alignment_a, genomes_a)

        # store ortholog name retrieved from filename
        ortholog = os.path.basename(sico_file).split('.')[0]
        instance.values[ORTHOLOG] = ortholog

        # add codeml_values to clade_calcs instance values
        instance.values.update(codeml_values)

        # add phipack values for this file
        instance.values.update(phipack_values[sico_file])

        # add COG digits and letters
        _extract_cog_digits_and_letters(instance)

        # add SFS related values
        _codon_site_freq_spec(instance)

        # add additional deduced calculation
        _add_combined_calculations(instance)

        # store the clade_calc values
        calculations.append(instance)

    # calculcate mean and averages
    max_nton = len(genome_ids_a) // 2
    sum_stats, mean_stats = _calculcate_mean_and_averages(calculations, max_nton)

    # neutrality index calculation and bootstrapping
    ni_stats, ni_lower_stats, ni_upper_stats = _neutrality_indices(calculations)

    # finally append statistics to calculations so they show up in file
    calculations.extend((sum_stats, mean_stats, ni_stats, ni_lower_stats, ni_upper_stats))

    return calculations


def run_calculations(genomes_a_file,
                     genomes_b_file,
                     sico_files,
                     table_a_dest,
                     table_b_dest):
    '''Perform all calculations as requested through command line arguments'''
    # parse genomes in genomes_x_files
    genome_ids_a, common_prefix_a = _extract_genome_ids_and_common_prefix(genomes_a_file)
    genome_ids_b, common_prefix_b = _extract_genome_ids_and_common_prefix(genomes_b_file)

    # calculate phipack values for combined aligments once
    if DEBUG:
        # PhiPack is SLOW, so when debugging just return zero
        phipack_values = {sico_file:
                          defaultdict(int)
                          for sico_file in sico_files}
    else:
        phipack_dir = tempfile.mkdtemp(prefix='phipack_')
        phipack_values = {sico_file:
                          run_phipack(phipack_dir, sico_file)
                          for sico_file in sico_files}
        shutil.rmtree(phipack_dir)

    # per table calculations
    if 1 < len(genome_ids_a):
        calculations_ab = _table_calculations(genome_ids_a, genome_ids_b, sico_files, phipack_values)
        _write_to_file(table_a_dest,
                       genome_ids_a, genome_ids_b,
                       common_prefix_a, common_prefix_b,
                       calculations_ab)
    else:
        with open(table_a_dest, 'w') as write_handle:
            write_handle.write('#At least two genomes are needed to calculate diversity, not ' + str(len(genome_ids_a)))

    if 1 < len(genome_ids_b):
        calculations_ba = _table_calculations(genome_ids_b, genome_ids_a, sico_files, phipack_values)
        _write_to_file(table_b_dest,
                       genome_ids_b, genome_ids_a,
                       common_prefix_b, common_prefix_a,
                       calculations_ba)
    else:
        with open(table_b_dest, 'w') as write_handle:
            write_handle.write('#At least two genomes are needed to calculate diversity, not ' + str(len(genome_ids_b)))


def _every_other_codon_alignments(alignment):
    """Separate alignment into separate alignments per codon, to get independent axis when graphing data."""
    # Calculate sequence_length to use when splitting MSA into codons
    sequence_lengths = len(alignment[0]) - len(alignment[0]) % 3
    alignment_codons = [alignment[:, index:index + 3] for index in range(0, sequence_lengths, 3)]

    def _concat_codons_to_alignment(codons):
        """Concatenate codons as Bio.Align.MultipleSeqAlignment to one another to create a composed MSA."""
        alignment = codons[0]
        for codon in codons[1:]:
            alignment += codon
        return alignment

    # Odd alignments are the sum of the odd codons
    ali_odd = _concat_codons_to_alignment([codon for index, codon in enumerate(alignment_codons, 1) if index % 2 == 1])
    # Even alignments are the sum of the even codons
    ali_even = _concat_codons_to_alignment([codon for index, codon in enumerate(alignment_codons, 1) if index % 2 == 0])
    return ali_odd, ali_even


def _split_by_odd_even_codons(sico_files):
    '''Split each sequence in each sico file by odd and even codons and write them out to separate files.'''
    odd_sico_files = []
    even_sico_files = []

    for sico_file in sico_files:
        # split alignments and store them in separate files
        alignment = AlignIO.read(sico_file, 'fasta')
        odd_alignment, even_alignment = _every_other_codon_alignments(alignment)

        # compose new unique filenames
        dirname, basename = os.path.split(sico_file)
        odd_file = os.path.join(dirname, 'odd_' + basename)
        even_file = os.path.join(dirname, 'even_' + basename)

        # store the split alignments in separate files
        AlignIO.write(odd_alignment, odd_file, 'fasta')
        AlignIO.write(even_alignment, even_file, 'fasta')

        # safe the references to the new files
        odd_sico_files.append(odd_file)
        even_sico_files.append(even_file)

    return odd_sico_files, even_sico_files


def _prepare_calculations(genomes_a_file,
                          genomes_b_file,
                          sicozip_file,
                          table_a_dest,
                          table_b_dest,
                          append_odd_even=False):
    '''Unzip sico_files, and if needed create temporary files for the odd/even only codons.'''
    if append_odd_even:
        # prepend file makeup when odd/even table are also added
        _write_intro_to_file(table_a_dest)
        _write_intro_to_file(table_b_dest)

    # extract ortholog files from sicozip
    rundir = tempfile.mkdtemp(prefix='calculations_')
    sico_files = extract_archive_of_files(sicozip_file, create_directory('sicos', inside_dir=rundir))

    # perform normal calculation
    run_calculations(genomes_a_file, genomes_b_file, sico_files, table_a_dest, table_b_dest)

    # separate calculations for odd and even tables
    if append_odd_even:
        odd_sico_files, even_sico_files = _split_by_odd_even_codons(sico_files)
        run_calculations(genomes_a_file, genomes_b_file, odd_sico_files, table_a_dest, table_b_dest)
        run_calculations(genomes_a_file, genomes_b_file, even_sico_files, table_a_dest, table_b_dest)

    # clean up
    shutil.rmtree(rundir)

def main(argv=None):  # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv[1:]

    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __doc__.split("\n")[1]
    program_usage = '''%s

  Created by Tim te Beek on %s.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_usage, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %(default)s]")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)

        def test_file_readable(argument, mode=os.R_OK):
            if os.access(argument, mode):
                return argument
            raise ArgumentTypeError('File {} is not {}'.format(argument, mode))

        # Arguments specific to calculations
        parser.add_argument('--genomes-a', nargs=1, type=test_file_readable, required=True,
                            help='Tab separated values file with Genome IDs of clade A')
        parser.add_argument('--genomes-b', nargs=1, type=test_file_readable, required=True,
                            help='Tab separated values file with Genome IDs of clade B')
        parser.add_argument('--sico-zip', nargs=1, type=test_file_readable, required=True,
                            help='Zip archive containing Single Copy Ortholog files')
        parser.add_argument('--table-a', nargs=1, default='table-a.tsv',
                            help='Destination output file path for comparison of clade A with clade B')
        parser.add_argument('--table-b', nargs=1, default='table-b.tsv',
                            help='Destination output file path for comparison of clade B with clade A')

        parser.add_argument('-a', '--append-odd-even', action='store_true',
                            help='append separate tables calculated for odd and even codons of ortholog alignments (default: False)')

        # Process arguments
        args = parser.parse_args(argv)

        if args.verbose > 0:
            print("Verbose mode on")
            logging.root.setLevel(logging.DEBUG)

        # perform the calculations
        _prepare_calculations(args.genomes_a[0],
                              args.genomes_b[0],
                              args.sico_zip[0],
                              args.table_a[0],
                              args.table_b[0],
                              args.append_odd_even)

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0

if __name__ == "__main__":
    sys.exit(main())
