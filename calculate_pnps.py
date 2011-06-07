#!/usr/bin/env python
"""Module to calculate pn ps."""
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable
from divergence import parse_options
import logging as log
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

def _perform_calculations(alignment):
    """Perform actual calculations on the alignment to determine pN, pS, SFS & the number of ignored cases per SICO."""
    #Print codons for easy debugging
    for seqr in alignment:
        seq = str(seqr.seq)
        log.info('  '.join(seq[idx:idx + 3] for idx in range(0, len(seq), 3)))

    synonymous_polymorphisms = 0
    synonymous_sfs = {}
    non_synonymous_polymorphisms = 0
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
        lsfs = dict((ntimes, psu_values.count(ntimes)) for ntimes in set(psu_values))
        #Deduct one for the reference_allele_count, which should not count towards the SFS
        lsfs[reference_allele_count] = lsfs[reference_allele_count] - 1
        if lsfs[reference_allele_count] == 0:
            del lsfs[reference_allele_count]

        if synonymous:
            #If all polymorphisms encode for the same AA, we have multiple synonymous polymorphisms, where:
            #2 nucleotides = 1 polymorphism, 3 nucleotides = 2 polymorphisms, 4 nucleotides = 3 polymorphisms
            synonymous_polymorphisms += sum(lsfs.values())

            #Debug log statement
            log.info('local SFS {0}'.format(lsfs))

            #Increment count of minor allele occupations in dictionary that tracks this per alignment
            for maf, count in lsfs.iteritems():
                prev_occupations = synonymous_sfs.get(maf, 0)
                synonymous_sfs[maf] = prev_occupations + count
            continue

        if not synonymous:
            if len(polymorph_site_usage) == len(translation_usage):
                #If all polymorphisms encode for different AA, we have multiple non-synonymous polymorphisms, where:
                #2 nucleotides = 1 polymorphism, 3 nucleotides = 2 polymorphisms, 4 nucleotides = 3 polymorphisms
                non_synonymous_polymorphisms += sum(lsfs.values())

                #Debug log statement
                log.info('local SFS {0}'.format(lsfs))

                #Increment count of minor allele occupations in dictionary that tracks this per alignment
                for maf, count in lsfs.iteritems():
                    prev_occupations = non_synonymous_sfs.get(maf, 0)
                    non_synonymous_sfs[maf] = prev_occupations + count
                continue
            else:
                #Some, but not all polymorphisms encode for different AA, making it unclear how this should be scored
                mixed_synonymous_polymorphisms += 1
                continue

    log.info('synonymous_polymorphisms %i', synonymous_polymorphisms)
    log.info('synonymous_sfs %s', str(synonymous_sfs))
    log.info('non_synonymous_polymorphisms %i', non_synonymous_polymorphisms)
    log.info('non_synonymous_sfs %s', str(non_synonymous_sfs))
    log.info('mixed_synonymous_polymorphisms %i', mixed_synonymous_polymorphisms)
    log.info('multiple_site_polymorphisms %i', multiple_site_polymorphisms)

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
