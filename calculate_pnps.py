#!/usr/bin/env python
"""Module to calculate pn ps."""
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable
from divergence import parse_options
from itertools import product
import logging as log
import re
import sys

def calculate_pnps(genome_ids_a, genome_ids_b, sico_files):
    """Calculate pN, pS, pN/pS & SFS values for all sico_files per clade, and write out statistics per SICO."""
    #For each alignment create separate alignments for clade A & clade B genomes 
    alignments = (AlignIO.read(sico_file, 'fasta') for sico_file in sico_files)
    for ali in alignments:
        alignment_a = MultipleSeqAlignment(seqr for seqr in ali if seqr.id.split('|')[0] in genome_ids_a)
        alignment_b = MultipleSeqAlignment(seqr for seqr in ali if seqr.id.split('|')[0] in genome_ids_b)

        #Calculate pn ps for the subaligments of each clade
        _perform_calculations(alignment_a)
        _perform_calculations(alignment_b)

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

def _perform_calculations(alignment):
    """Perform actual calculations on the alignment to determine pN, pS, SFS & the number of ignored cases per SICO."""
    #Print codons for easy debugging
    for seqr in alignment:
        seq = str(seqr.seq)
        log.info('  '.join(seq[idx:idx + 3] for idx in range(0, len(seq), 3)))

    synonymous_sfs = {}
    four_fold_syn_sfs = {}
    non_synonymous_sfs = {}
    mixed_synonymous_polymorphisms = 0
    multiple_site_polymorphisms = 0

    #Calculate range_end here so we can handle alignments that are not multiples of three
    range_end = len(alignment[0]) - len(alignment[0]) % 3
    #Split into codon_alignments
    codon_alignments = [alignment[:, index:index + 3] for index in range(0, range_end, 3)]
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

        #Mutations are synonymous when all codons encode the same AA, and there are no gaps 
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
            continue

        #Determine if only one site is polymorphic by using boolean xor and not all
        single_site_polymorphism = site1_polymorphic ^ site2_polymorphic ^ site3_polymorphic and not all(polymorphisms)

        #Debug print statement        
        log.info('1:{0}  2:{1}  3:{2}  Translations:{3}\tSynonymous:{4}\tSingle-site:{5}' \
        .format(site1_usage, site2_usage, site3_usage, translation_usage, synonymous, single_site_polymorphism))

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
        local_sfs = dict((ntimes, psu_values.count(ntimes)) for ntimes in set(psu_values))
        #Deduct one for the reference_allele_count, which should not count towards the SFS
        local_sfs[reference_allele_count] = local_sfs[reference_allele_count] - 1
        if local_sfs[reference_allele_count] == 0:
            del local_sfs[reference_allele_count]

        if synonymous:
            #If all polymorphisms encode for the same AA, we have multiple synonymous polymorphisms, where:
            #2 nucleotides = 1 polymorphism, 3 nucleotides = 2 polymorphisms, 4 nucleotides = 3 polymorphisms

            #Debug log statement
            log.info('local SFS: {0}'.format(local_sfs))

            #Update synonymous SFS by adding values from local SFS
            for maf, count in local_sfs.iteritems():
                prev_occupations = synonymous_sfs.get(maf, 0)
                synonymous_sfs[maf] = prev_occupations + count

            #Codon is four fold degenerate if it matches any pattern in FOUR_FOLD_DEGENERATE_PATTERNS
            if site3_polymorphic:
                codon = site1[0] + site2[0] + site3[0]
                for pattern in FOUR_FOLD_DEGENERATE_PATTERNS:
                    if re.match(pattern, codon):
                        #Update four fold degenerate SFS by adding values from local SFS
                        for maf, count in local_sfs.iteritems():
                            prev_occupations = four_fold_syn_sfs.get(maf, 0)
                            four_fold_syn_sfs[maf] = prev_occupations + count
            continue

        if not synonymous:
            if len(polymorph_site_usage) == len(translation_usage):
                #If all polymorphisms encode for different AA, we have multiple non-synonymous polymorphisms, where:
                #2 nucleotides = 1 polymorphism, 3 nucleotides = 2 polymorphisms, 4 nucleotides = 3 polymorphisms

                #Debug log statement
                log.info('local SFS: {0}'.format(local_sfs))

                #Update non synonymous SFS by adding values from local SFS
                for maf, count in local_sfs.iteritems():
                    prev_occupations = non_synonymous_sfs.get(maf, 0)
                    non_synonymous_sfs[maf] = prev_occupations + count
                continue
            else:
                #Some, but not all polymorphisms encode for different AA, making it unclear how this should be scored
                mixed_synonymous_polymorphisms += 1
                continue

    #Use Site Frequency Spectrum to calculate the number of synonymous and non synonymous polymorphisms
    #Note: this requires a complete SFS across synonymous & non_synonymous polymorphisms, be careful when updating code
    synonymous_polymorphisms = sum(synonymous_sfs.values())
    four_fold_syn_polymorphisms = sum(four_fold_syn_sfs.values())
    non_synonymous_polymorphisms = sum(non_synonymous_sfs.values())

    log.info('\n')

    log.info('synonymous_polymorphisms %i', synonymous_polymorphisms)
    log.info('synonymous_sfs:\t%s\n', _sfs2str(synonymous_sfs))

    log.info('four fold synonymous_polymorphisms %i', four_fold_syn_polymorphisms)
    log.info('four fold synonymous_sfs:\t%s\n', _sfs2str(four_fold_syn_sfs))

    log.info('non_synonymous_polymorphisms %i', non_synonymous_polymorphisms)
    log.info('non_synonymous_sfs:\t%s\n', _sfs2str(non_synonymous_sfs))

    log.info('mixed_synonymous_polymorphisms %i', mixed_synonymous_polymorphisms)
    log.info('multiple_site_polymorphisms %i', multiple_site_polymorphisms)

    #1. gene name
    gene_name = 'xyz'

    #2. number of strains (this will typically be the same for all genes)
    number_of_strains = len(alignment)

    #3. number of non-synonymous sites (see below for calculation) = LnP
    #TODOnumber_of_non_syn_sites = ?

    #4. SFS for non-synonymous polymorphsims
    non_synonymous_sfs

    #5. The sum of SFS for non-synonymous polymorphisms = Pn
    pn = sum(non_synonymous_sfs.values())

    #6. number of synonymous sites (see below for calculation) = LsP
    #TODO number_of_non_syn_sites = ?

    #7. SFS for synonymous prolymorphisms
    synonymous_sfs

    #8. The sum of the SFS for synonymous polymorphisms = Ps
    ps = sum(synonymous_sfs.values())

    #9. number of 4-fold synonymous sites = L4


    #10. SFS for 4-fold synonymous polymorphisms


    #11. Sum of the SFS for 4-folds = P4


    #12. number of non-synonymous sites for divergence from PAML (this might be different to 3, because this will come from two randomly chosen sequences) = LnD


    #13. number of non-synonymous substitutions from PAML = Dn


    #14. number of synonymous sites from PAML = LsD


    #15. number of synonymous substitutions from PAML = Ds


    #16. Direction of selection = Dn/(Dn+Ds) - Pn/(Pn+Ps)


    #I would calculate LnP and LsP as follows:
    #
    #LnP = L * LnD/(LnD+LsD)
    #LsP = L * LsD/(LnD+LsD)
    #
    #where L is the length of the sequence from which the polymorphism data is taken (i.e. 3* no of codons used)
    #
    #Finally, we might a line which gives the sum for columns 3 to 15 plus the mean of column 16

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: calculate_pnps.py 
--alignment=FILE    alignment file in FASTA format
"""
    options = ['alignment']
    (alignment_file,) = parse_options(usage, options, args)

    alignment = AlignIO.read(alignment_file, 'fasta')
    _perform_calculations(alignment)

if __name__ == '__main__':
    main(sys.argv[1:])
