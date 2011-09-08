'''
Created on Jul 7, 2011

@author: tbeek
'''

from __future__ import division
from Bio import SeqIO
from Bio.Data import CodonTable
from divergence import CODON_TABLE_ID
from divergence.select_taxa import _parse_genomes_table, download_genome_files
from operator import itemgetter
from random import choice
import logging
import sys

def _sample_genomes(sample_size = 100):
    """Take random genomes from the microbial genomes table, without duplicating the first part of the name."""
    #We want random unique species, so store the previous organism names
    previous_names = set()

    #Parse complete microbial genomes table to produce dictionaries containing column headers mapped to values
    genomes = _parse_genomes_table(require_refseq = True)

    while len(previous_names) < sample_size:
        genome = choice(genomes)

        if genome['Super Kingdom'] != 'Bacteria':#Archaea / Bacteria
            continue

        #Skip some known troublesome cases
        if genome['RefSeq project ID'] in ('57665', '57721', '57727', '57755', '57757', '57853', '58069', '58599', \
                                           '58901', '58941', '59133', '62947', '68249', '59189', '58407', '58709'):
            continue

        firstname = genome['Organism Name'].split()[0]

        #Skip genome if organism firstname was previously selected
        if firstname in previous_names:
            continue

        #Yield organism
        logging.info('Selected:\t%s %s', genome['RefSeq project ID'], genome['Organism Name'])
        previous_names.add(firstname)
        yield genome

def _read_genbank_file(genbank_file):
    """Parse GenBank file using BioPython. Future extension point for sanity checks on input."""
    print genbank_file
    return SeqIO.read(genbank_file, 'genbank')

#Using the standard NCBI Bacterial, Archaeal and Plant Plastid Code translation table (11).
BACTERIAL_CODON_TABLE = CodonTable.unambiguous_dna_by_id.get(CODON_TABLE_ID)

def _extract_coding_sequences(gbk_record):
    """Extract coding sequences from genbank file from positive strand and return start, end and stopcodon."""
    #Skip if the first feature is not source
    first_feature = gbk_record.features[0]
    if first_feature.type != 'source':
        print >> sys.stderr, 'Expected first feature to be source, but was:', first_feature
        return

    #Skip if the source feature tells us it's from a plasmid 
    if 'plasmid' in first_feature.qualifiers:
        #Skipping plasmids for now
        raise StopIteration

    coding_features = (gb_feature
                        for gb_feature in gbk_record.features
                        #Skip any non coding sequence features or pseudo (non-functional version) CDS
                        if gb_feature.type == 'CDS' and not 'pseudo' in gb_feature.qualifiers)

    #Select only the coding sequences from all feature records: Bio.SeqFeature
    for cds_feature in coding_features:
        #Handle positive and negative strand separately
        if cds_feature.strand == 1:
            cds_start = cds_feature.location.start.position
            cds_end = cds_feature.location.end.position
        elif cds_feature.strand == -1:
            cds_start = 0 - cds_feature.location.end.position
            cds_end = 0 - cds_feature.location.start.position
        else:
            logging.warn('Expected strand value %i encountered in %s of %s', cds_feature.strand, cds_feature.id,
                         gbk_record.annotations["organism"])
            continue

        last_codon = str(cds_feature.extract(gbk_record.seq)[-3:])
        #Yield start position, end position and last codon if stopcodon
        if last_codon in BACTERIAL_CODON_TABLE.stop_codons:
            yield cds_start, cds_end, last_codon

def _group_cds_by_operons(coding_sequences, operon_distance = 300):
    """Group coding sequences into operons based on how close the start of a gene is to the end of the previous gene."""
    #Recurring method to retrieve end of coding sequence feature
    cds_start_func = itemgetter(0)
    cds_end_func = itemgetter(1)

    operon = []
    #This sorted call retrieves all coding sequences in one go
    #This code does not yet take into account the positive and negative strand handling introduced in _extract_cds above
    for cds in sorted(coding_sequences, key = cds_end_func):
        #Append coding sequence to operon if it's start is smaller than the current operon's largest end + distance 
        if not operon or cds_start_func(cds) < max(cds_end_func(cds2) for cds2 in operon) + operon_distance:
            operon.append(cds)
        else:
            #Return operons larger than one
            if 1 < len(operon):
                yield operon
            #Start new operon with current cds as first
            operon = [cds]

def _usage_per_operon_position(operons):
    """Calculate usage for each of the positions"""
    stopcodon_usage_per_operon_position = {}
    largest_operon = max(len(operon) for operon in operons)
    for operon_position in range(0, largest_operon):
        position_stopcodons = [operon[operon_position][2]
                               for operon in operons if operon_position < len(operon)]
        position_stopcodon_counts = dict((codon, position_stopcodons.count(codon))
                                         for codon in set(position_stopcodons))
        stopcodon_usage_per_operon_position[operon_position] = position_stopcodon_counts
    return stopcodon_usage_per_operon_position

def __main__():
    """Main method called when run from terminal."""
    #Select random genomes
    genomes = _sample_genomes(100)

    #Download genbank files for genomes in the background
    from multiprocessing import Pool
    pool = Pool()
    ft_genome_files = [(genome['Organism Name'], pool.apply_async(download_genome_files, (genome,)))
                       for genome in genomes]
    genome_files = ((name, (pairs[1]
                            for pairs in future.get()))
                    for name, future in ft_genome_files)

    #Extract coding sequences from genome files
    genome_gbk_records = ((name, (_read_genbank_file(genbank_file)
                                  for genbank_file in genbank_files))
                          for name, genbank_files in genome_files)
    genome_stopcodons = ((name, [cds_tuple[2]
                                 for gbk_record in gbk_records
                                 for cds_tuple in _extract_coding_sequences(gbk_record)])
                         for name, gbk_records in genome_gbk_records)

    #Stopcodons per genome
    usages = ((name, dict((codon, stopcodons.count(codon))
                   for codon in set(stopcodons)))
              for name, stopcodons in genome_stopcodons
              if stopcodons)

    #Write out usages
    with open('stopcodon-usage-across-species.tsv', mode = 'w') as write_handle:
        write_handle.write('name\tTAA\t%\tTAG\t%\tTGA\t%\t\n')
        for name, usage in usages:
            total = sum(usage.values())
            taa = usage.get('TAA', 0)
            tag = usage.get('TAG', 0)
            tga = usage.get('TGA', 0)
            write_handle.write('\t'.join(str(item) for item in (name,
                                                                taa, taa / total * 100,
                                                                tag, tag / total * 100,
                                                                tga, tga / total * 100)) + '\n')

    print 'Done!'

if __name__ == '__main__':
    __main__()
