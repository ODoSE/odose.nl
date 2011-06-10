#!/usr/bin/env python
"""Module to calculate pn ps."""

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable
from divergence import parse_options, create_directory
from divergence.run_codeml import run_codeml, parse_codeml_output
from itertools import product
import logging as log
import os.path
import re
import sys
import tempfile

def calculate_pnps(genome_ids_a, genome_ids_b, sico_files):
    """Calculate pN, pS, pN/pS & SFS values for all sico_files per clade, and write out statistics per SICO."""
    run_dir = tempfile.mkdtemp(prefix = 'calculate_pnps')

    #For each alignment create separate alignments for clade A & clade B genomes
    calculations_a_file = os.path.join(run_dir, 'calculations_a.tsv')
    calculations_b_file = os.path.join(run_dir, 'calculations_b.tsv')
    _write_output_file_header(calculations_a_file)
    _write_output_file_header(calculations_b_file)

    for sico_file in sico_files:
        ali = AlignIO.read(sico_file, 'fasta')
        alignment_a = MultipleSeqAlignment(seqr for seqr in ali if seqr.id.split('|')[0] in genome_ids_a)
        alignment_b = MultipleSeqAlignment(seqr for seqr in ali if seqr.id.split('|')[0] in genome_ids_b)

        #1. gene name
        basename = os.path.split(sico_file)[1].split('.')[0]

        #Run codeml to calculate values for dn & ds
        codeml_file = run_codeml(create_directory(basename, inside_dir = run_dir), alignment_a, alignment_b)
        value_dict = parse_codeml_output(codeml_file)

        #Perform calculations for subaligments of each clade, if clade has more than one sequence; skipping outliers
        if len(alignment_a):
            comp_values = _perform_calculations(alignment_a, value_dict)
            _append_statistics(calculations_a_file, basename, comp_values)
        if len(alignment_b):
            comp_values = _perform_calculations(alignment_b, value_dict)
            _append_statistics(calculations_b_file, basename, comp_values)

    #Finally, we might need a line which gives the sum for columns 3 to 15 plus the mean of column 16

    #Clean up
    #TODO After moving files to proper paths outside removed directory: shutil.rmtree(run_dir)

#Using the standard NCBI Bacterial, Archaeal and Plant Plastid Code translation table (11).
BACTERIAL_CODON_TABLE = CodonTable.unambiguous_dna_by_id.get(11)

def _four_fold_degenerate_patterns():
    """Find patterns of 4-fold degenerate codons, wherein all third site substitutions code for the same amino acid."""
    letters = BACTERIAL_CODON_TABLE.nucleotide_alphabet.letters
    #Any combination of letters of length two
    for site12 in [''.join(prod) for prod in product(letters, repeat = 2)]:
        #4-fold when the length of the unique encoded amino acids for all possible third site nucleotides is exactly 1 
        if 1 == len(set([BACTERIAL_CODON_TABLE.forward_table.get(site12 + site3) for site3 in letters])):
            #Add regular expression pattern to the set of patterns
            yield '{0}[{1}]'.format(site12, letters)

FOUR_FOLD_DEGENERATE_PATTERNS = set(_four_fold_degenerate_patterns())

def _sfs2str(sfs):
    """Convert Site Frequency Spectrum to human readable textual representation."""
    def _generator(sfs):
        """Return successive values for - and names of multipleton(s)."""
        multiple_names = {1:'singleton', 2:'doubleton', 3:'tripleton', 4:'quadrupleton', 5:'quintupleton'}
        for xton in sorted(sfs.keys()):
            name = multiple_names.get(xton, '"{0}"-ton'.format(xton))
            yield '{0} {1}'.format(sfs[xton], name) + ('s' if 1 < sfs[1] else '')
    return ', '.join(_generator(sfs))

def _perform_calculations(alignment, codeml_values):
    """Perform actual calculations on the alignment to determine pN, pS, SFS & the number of ignored cases per SICO."""
    synonymous_sfs = {}
    four_fold_syn_sfs = {}
    non_synonymous_sfs = {}
    four_fold_synonymous_sites = 0
    mixed_synonymous_polymorphisms = 0
    multiple_site_polymorphisms = 0

    #Calculate alignment_length here so we can handle alignments that are not multiples of three
    alignment_length = len(alignment[0]) - len(alignment[0]) % 3
    #Split into codon_alignments
    codon_alignments = [alignment[:, index:index + 3] for index in range(0, alignment_length, 3)]
    for codon_alignment in codon_alignments:
        #Get string representations of codons for simplicity 
        codons = [str(seqr.seq) for seqr in codon_alignment]

        #As per AEW: ignore codons with gaps, and codons with unresolved bases: Basically anything but ACGT
        if 0 < len(''.join(codons).translate(None, 'ACGTactg')):
            continue

        #Stop codons should have already been removed, but lets be sure anyway
        codons_nonstop = [codon for codon in codons if codon not in BACTERIAL_CODON_TABLE.stop_codons]

        #Retrieve translations of codons now that inconclusive & stop-codons have been removed
        translations = [BACTERIAL_CODON_TABLE.forward_table.get(codon) for codon in codons_nonstop]

        #Count unique translations across strains
        translation_usage = dict((aa, translations.count(aa)) for aa in set(translations))

        #Mutations are synonymous when all codons encode the same AA, and there are no skipped codons
        synonymous = len(translation_usage) == 1 and len(translations) == len(codon_alignment)

        #Retrieve nucleotides per site within the codon 
        site1 = [nucl for nucl in codon_alignment[:, 0]]
        site2 = [nucl for nucl in codon_alignment[:, 1]]
        site3 = [nucl for nucl in codon_alignment[:, 2]]

        #Count occurrences of distinct nucleotides across strains
        site1_usage = dict((nucl, site1.count(nucl)) for nucl in set(site1))
        site2_usage = dict((nucl, site2.count(nucl)) for nucl in set(site2))
        site3_usage = dict((nucl, site3.count(nucl)) for nucl in set(site3))

        #Sites are polymorphic if they contain more than one nucleotide
        site1_polymorphic = 1 < len(site1_usage)
        site2_polymorphic = 1 < len(site2_usage)
        site3_polymorphic = 1 < len(site3_usage)
        polymorphisms = site1_polymorphic, site2_polymorphic, site3_polymorphic

        #Continue with next codon if none of the sites is polymorphic
        if not any(polymorphisms):
            #But do increase the number of 4-fold synonymous sites if the pattern matches
            codon = codons_nonstop[0]
            for pattern in FOUR_FOLD_DEGENERATE_PATTERNS:
                if re.match(pattern, codon):
                    #Increase by one, as this site is for fold degenerate, even if it is not polymorphic
                    four_fold_synonymous_sites += 1
            continue

        #Determine if only one site is polymorphic by using boolean xor and not all
        single_site_polymorphism = site1_polymorphic ^ site2_polymorphic ^ site3_polymorphic and not all(polymorphisms)

        #Skip multiple site polymorphisms, but do keep a count of how many we encounter
        if not single_site_polymorphism:
            multiple_site_polymorphisms += 1
            continue

        #Determine which site_usage is the single site polymorphism
        polymorph_site_usage = site1_usage if site1_polymorphic else site2_usage if site2_polymorphic else site3_usage

        #Find the 'reference' nucleotide as (one of) the most occurring occupations in this site, so we can -1 later
        psu_values = polymorph_site_usage.values()
        reference_allele_count = max(psu_values)

        #Calculate the local site frequency spectrum, to be added to the gene-wide SFS later
        #We'll be using Site Frequency Spectrum to calculate the number of synonymous and non synonymous polymorphisms
        #Note: this requires a complete SFS across synonymous & non_synonymous polymorphisms, be careful when updating
        local_sfs = dict((ntimes, psu_values.count(ntimes)) for ntimes in set(psu_values))
        #Deduct one for the reference_allele_count, which should not count towards the SFS
        local_sfs[reference_allele_count] = local_sfs[reference_allele_count] - 1
        #Remove empty value as possible result of the above decrement operation
        if local_sfs[reference_allele_count] == 0:
            del local_sfs[reference_allele_count]

        def _update_sfs_with_local_sfs(sfs, local_sfs):
            """Add values from local_sfs to gene-wide sfs"""
            for maf, count in local_sfs.iteritems():
                prev_occupations = sfs.get(maf, 0)
                sfs[maf] = prev_occupations + count

        if synonymous:
            #If all polymorphisms encode for the same AA, we have multiple synonymous polymorphisms, where:
            #2 nucleotides = 1 polymorphism, 3 nucleotides = 2 polymorphisms, 4 nucleotides = 3 polymorphisms

            #Update synonymous SFS by adding values from local SFS
            _update_sfs_with_local_sfs(synonymous_sfs, local_sfs)

            #Codon is four fold degenerate if it matches any pattern in FOUR_FOLD_DEGENERATE_PATTERNS
            if site3_polymorphic:
                codon = codons_nonstop[0]
                for pattern in FOUR_FOLD_DEGENERATE_PATTERNS:
                    if re.match(pattern, codon):
                        #Update four fold degenerate SFS by adding values from local SFS
                        _update_sfs_with_local_sfs(four_fold_syn_sfs, local_sfs)
                        #Increase the number of four_fold synonymous sites here as well
                        four_fold_synonymous_sites += 1
        else: #not synonymous
            if len(polymorph_site_usage) == len(translation_usage):
                #If all polymorphisms encode for different AA, we have multiple non-synonymous polymorphisms, where:
                #2 nucleotides = 1 polymorphism, 3 nucleotides = 2 polymorphisms, 4 nucleotides = 3 polymorphisms

                #Update non synonymous SFS by adding values from local SFS
                _update_sfs_with_local_sfs(non_synonymous_sfs, local_sfs)
            else:
                #Some, but not all polymorphisms encode for different AA, making it unclear how this should be scored
                mixed_synonymous_polymorphisms += 1

    log.info('%i Polymorphisms with mixed synonymous and non-synonymous polymorphisms', mixed_synonymous_polymorphisms)
    log.info('%i of %i codons have polymorphisms in multiple sites', multiple_site_polymorphisms, alignment_length / 3)

    #Compute combined values from the above counted statistics 
    return _compute_values_from_statistics(alignment, alignment_length, codeml_values,
                                           synonymous_sfs, non_synonymous_sfs,
                                           four_fold_syn_sfs, four_fold_synonymous_sites)

def _compute_values_from_statistics(alignment, alignment_length, codeml_values, synonymous_sfs, non_synonymous_sfs,
                                    four_fold_syn_sfs, four_fold_synonymous_sites):
    """Compute values (mostly from the site frequency spectra) that we'll output according to provided formula's."""
    #12. number of non-synonymous sites for divergence from PAML
    #(this might be different to 3, because this will come from two randomly chosen sequences) = LnD
    paml_non_synonymous_sites = codeml_values['N']

    #14. number of synonymous sites from PAML = LsD
    paml_synonymous_sites = codeml_values['S']

    #Dictionary to hold computed values
    calc_values = {}

    #Add all values in codeml_values to calc_values
    calc_values.update(codeml_values)

    #2. number of strains (this will typically be the same for all genes)
    calc_values['strains'] = len(alignment)

    #3. number of non-synonymous sites (see below for calculation) = LnP = L * LnD/(LnD+LsD)
    #where L is the length of the sequence from which the polymorphism data is taken (i.e. 3* no of codons used)
    calc_values['non-synonymous sites'] = alignment_length * paml_non_synonymous_sites / (paml_non_synonymous_sites +
                                                                                          paml_synonymous_sites)

    #4. SFS for non-synonymous polymorphsims
    calc_values['non-synonymous sfs'] = non_synonymous_sfs

    #5. The sum of SFS for non-synonymous polymorphisms = Pn
    calc_values['non-synonymous polymorphisms'] = non_synonymous_polymorphisms = sum(non_synonymous_sfs.values())

    #6. number of synonymous sites (see below for calculation) = LsP = L * LsD/(LnD+LsD)
    #where L is the length of the sequence from which the polymorphism data is taken (i.e. 3* no of codons used)
    calc_values['synonymous sites'] = alignment_length * paml_synonymous_sites / (paml_non_synonymous_sites +
                                                                                  paml_synonymous_sites)

    #7. SFS for synonymous polymorphisms
    calc_values['synonymous sfs'] = synonymous_sfs

    #8. The sum of the SFS for synonymous polymorphisms = Ps
    calc_values['synonymous polymorphisms'] = synonymous_polymorphisms = sum(synonymous_sfs.values())

    #9. number of 4-fold synonymous sites = L4  (all degenerate sites, _including_ those that are not polymorphic)
    calc_values['4-fold synonymous sites'] = four_fold_synonymous_sites

    #10. SFS for 4-fold synonymous polymorphisms
    calc_values['4-fold synonymous sfs'] = four_fold_syn_sfs

    #11. Sum of the SFS for 4-folds = P4
    calc_values['4-fold synonymous polymorphisms'] = sum(four_fold_syn_sfs.values())

    #13. number of non-synonymous substitutions from PAML = Dn
    paml_non_synonymous_substitut = codeml_values['Dn']

    #15. number of synonymous substitutions from PAML = Ds
    paml_synonymous_substitutions = codeml_values['Ds']

    #16. Direction of selection = Dn/(Dn+Ds) - Pn/(Pn+Ps)
    calc_values['direction of selection'] = (paml_non_synonymous_substitut / (paml_non_synonymous_substitut +
                                                                                 paml_synonymous_substitutions)
                                            - non_synonymous_polymorphisms / (non_synonymous_polymorphisms +
                                                                              synonymous_polymorphisms))
    return calc_values

def _write_output_file_header(calculations_file):
    """Write header line for combined statistics file."""
    with open(calculations_file, mode = 'a') as append_handle:
        append_handle.write('{0}\t{1}\t'.format('orthologname',
                                                'strains'))
        append_handle.write('{0}\t{1}\t{2}\t'.format('non-synonymous sites',
                                                     'non-synonymous sfs',
                                                     'non-synonymous polymorphisms'))
        append_handle.write('{0}\t{1}\t{2}\t'.format('synonymous sites',
                                                     'synonymous sfs',
                                                     'synonymous polymorphisms'))
        append_handle.write('{0}\t{1}\t{2}\t'.format('4-fold synonymous sites',
                                                     '4-fold synonymous sfs',
                                                     '4-fold synonymous polymorphisms'))
        append_handle.write('{0}\t{1}\t{2}\t{3}\t'.format('N',
                                                          'Dn',
                                                          'S',
                                                          'Ds'))
        append_handle.write('{0}\n'.format('direction of selection'))


def _append_statistics(calculations_file, orthologname, comp_values):
    """Append statistics for individual """
    with open(calculations_file, mode = 'a') as append_handle:
        append_handle.write('{0}\t{1}\t'.format(orthologname,
                                                comp_values['strains']))
        append_handle.write('{0}\t{1}\t{2}\t'.format(comp_values['non-synonymous sites'],
                                                     comp_values['non-synonymous sfs'],
                                                     comp_values['non-synonymous polymorphisms']))
        append_handle.write('{0}\t{1}\t{2}\t'.format(comp_values['synonymous sites'],
                                                     comp_values['synonymous sfs'],
                                                     comp_values['synonymous polymorphisms']))
        append_handle.write('{0}\t{1}\t{2}\t'.format(comp_values['4-fold synonymous sites'],
                                                     comp_values['4-fold synonymous sfs'],
                                                     comp_values['4-fold synonymous polymorphisms']))
        append_handle.write('{0}\t{1}\t{2}\t{3}\t'.format(comp_values['N'],
                                                          comp_values['Dn'],
                                                          comp_values['S'],
                                                          comp_values['Ds']))
        append_handle.write('{0}\n'.format(comp_values['direction of selection']))

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: calculate_pnps.py 
--alignment=FILE    alignment file in FASTA format
"""
    options = ['alignment']
    (alignment_file,) = parse_options(usage, options, args)

    alignment = AlignIO.read(alignment_file, 'fasta')
    _perform_calculations(alignment, {'Dn':-1, 'Ds':-1, 'N':-1, 'S':-1})

if __name__ == '__main__':
    main(sys.argv[1:])
