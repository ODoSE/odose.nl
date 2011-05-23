#!/usr/bin/env python
"""Module to calculate pn ps."""
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable
from divergence import parse_options
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
        print '  '.join(seq[idx:idx + 3] for idx in range(0, len(seq), 3))

    synonymous_polymorphisms = 0
    non_synonymous_polymorphisms = 0
    minor_allele_occupations = {}
    mixed_synonymous_polymorphisms = 0
    multiple_site_polymorphisms = 0

    #Split into codon_alignments
    codon_alignments = [alignment[:, index:index + 3] for index in range(0, len(alignment[0]), 3)]
    for codon_alignment in codon_alignments:
        #Retrieve translations of codon, eliminating gap- & stopcodons
        codons = (str(seqr.seq) for seqr in codon_alignment)
        codons_nongap = (codon for codon in codons if '-' not in codon)
        #Stop codons should have already been removed, but lets be sure anyway
        codons_nonstop = (codon for codon in codons_nongap if codon not in BACTERIAL_CODON_TABLE.stop_codons)
        translations = [BACTERIAL_CODON_TABLE.forward_table.get(codon) for codon in codons_nonstop]

        #Count unique translations across strains
        translation_usage = dict((aa, translations.count(aa)) for aa in set(translations))

        #Mutations are synonymous when coded all codons encode the same AA 
        synonymous = len(translation_usage) == 1 and len(translations) == len(codon_alignment)

        #Retrieve nucleotides per site within the codon 
        site1 = [nucl for nucl in codon_alignment[:, 0]]
        site2 = [nucl for nucl in codon_alignment[:, 1]]
        site3 = [nucl for nucl in codon_alignment[:, 2]]

        #Count usages of distinct nucleotides across strains
        site1_usage = dict((nucl, site1.count(nucl)) for nucl in set(site1))
        site2_usage = dict((nucl, site2.count(nucl)) for nucl in set(site2))
        site3_usage = dict((nucl, site3.count(nucl)) for nucl in set(site3))

        #Sites are polymorphic if at least one site contains more than one nucleotide
        site1_polymorphic = 1 < len(site1_usage)
        site2_polymorphic = 1 < len(site2_usage)
        site3_polymorphic = 1 < len(site3_usage)
        polymorphisms = site1_polymorphic, site2_polymorphic, site3_polymorphic

        #Continue with next codon if none of the sites is polymorphic
        if not any(polymorphisms):
            continue

        #Determine if only one site is polymorphic by using boolean xor and all
        single_site_polymorphism = site1_polymorphic ^ site2_polymorphic ^ site3_polymorphic and not all(polymorphisms)

        #Debug print statement        
        print site1_usage, site2_usage, site3_usage, translation_usage, synonymous, single_site_polymorphism

        #Skip multiple site polymorphisms, but do keep a count of how many we encounter
        if not single_site_polymorphism:
            multiple_site_polymorphisms += 1
            continue

        #Determine which site_usage is the single site polymorphism
        polymorphism_usage = site1_usage if site1_polymorphic else site2_usage if site2_polymorphic else site3_usage

        #Determine in how many strains the minor allele occurs
        minor_allele_count = min(polymorphism_usage.values())

        if synonymous:
            #If all polymorphisms encode for the same AA, we have multiple synonymous polymorphisms, where:
            #2 nucleotides = 1 polymorphism, 3 nucleotides = 2 polymorphisms, 4 nucleotides = 3 polymorphisms
            synonymous_polymorphisms += len(polymorphism_usage) - 1

            #Increment count of minor allele occupations in dictionary that tracks this per alignment
            prev_occupations = minor_allele_occupations.get(minor_allele_count, 0)
            minor_allele_occupations[minor_allele_count] = prev_occupations + 1
            continue

        if not synonymous:
            if len(polymorphism_usage) == len(translation_usage):
                #If all polymorphisms encode for different AA, we have multiple non-synonymous polymorphisms, where:
                #2 nucleotides = 1 polymorphism, 3 nucleotides = 2 polymorphisms, 4 nucleotides = 3 polymorphisms
                non_synonymous_polymorphisms += len(polymorphism_usage) - 1

                #Increment count of minor allele occupations in dictionary that tracks this per alignment
                prev_occupations = minor_allele_occupations.get(minor_allele_count, 0)
                minor_allele_occupations[minor_allele_count] = prev_occupations + 1
                continue
            else:
                #Some, but not all polymorphisms encode for different AA, making it unclear how this should be scored
                mixed_synonymous_polymorphisms += 1
                continue

    #TODO Multiple polymorphisms at a single site add +2 or +3 to polymorphism counts, \ 
    #but minor_allele_count is only incremented once per site   

    print 'synonymous_polymorphisms', synonymous_polymorphisms
    print 'non_synonymous_polymorphisms', non_synonymous_polymorphisms
    print 'minor_allele_occupations', minor_allele_occupations
    print 'mixed_synonymous_polymorphisms', mixed_synonymous_polymorphisms
    print 'multiple_site_polymorphisms', multiple_site_polymorphisms

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: calculate_pnps.py 
--alignment=FILE    alignment file in FASTA format
"""
    options = ['alignment']
    alignment_file = parse_options(usage, options, args)

    alignment = AlignIO.read(alignment_file, 'fasta')
    _perform_calculations(alignment)

if __name__ == '__main__':
    main(sys.argv[1:])
