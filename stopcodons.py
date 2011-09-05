'''
Created on Jul 7, 2011

@author: tbeek

Hypothesis: There's a bias in the usage of stopcodons within operons. 
'''

from __future__ import division
from Bio import SeqIO
from Bio.Data import CodonTable
from divergence import CODON_TABLE_ID
from divergence.select_taxa import _parse_genomes_table, download_genome_files
from itertools import chain
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
                                           '58901', '58941', '59133', '62947', '68249', '59189', '58407'):
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
        position_stopcodons = [operon[operon_position][2] for operon in operons if operon_position < len(operon)]
        position_stopcodon_counts = dict((codon, position_stopcodons.count(codon)) for codon in set(position_stopcodons))
        stopcodon_usage_per_operon_position[operon_position] = position_stopcodon_counts
    return stopcodon_usage_per_operon_position

def __main__():
    """Main method called when run from terminal."""
    if True:
        #Select random genomes
        genomes = _sample_genomes(100)

        #Download genbank files for genomes in the background
        from multiprocessing import Pool
        pool = Pool()
        ft_genome_files = [pool.apply_async(download_genome_files, (genome,)) for genome in genomes]
        genome_files = (future.get() for future in ft_genome_files)
        genbank_files = (pairs[1] for files in genome_files for pairs in files)
    else:
        #Read current genbank from cache directory
        import glob
        #TODO Fix path
        genbank_files = list(glob.glob('/data/dev/workspace-python/lib-divergence/divergence-cache/refseq/*/*.gbk'))

    #Extract coding sequences from genome files
    gbk_records = (_read_genbank_file(genbank_file) for genbank_file in genbank_files)
    cds_per_genome = (_extract_coding_sequences(gbk_record) for gbk_record in gbk_records)

    #Group coding sequences from each genome into operons
    operons_per_genome = (tuple(_group_cds_by_operons(genome_cds)) for genome_cds in cds_per_genome)

    #Chain operons from all genomes
    all_operons = list(chain.from_iterable(operons_per_genome))

    #Per genome, calculate overall stopcodon usage
    with open('stopcodon-stats.tsv', mode = 'w', buffering = 1) as write_handle:
        #Print header
        header = '\t'.join(('Position',
                            'TAG',
                            'TAA',
                            'TGA'))
        print header
        write_handle.write(header + ' \n')

        #Calculate stopcodon usage
        stopcodon_usage_per_operon_position = _usage_per_operon_position(all_operons)

        #Write out usage of each of the three codons at each position of an operon
        for position in sorted(stopcodon_usage_per_operon_position.keys()):
            position_stopcodon_counts = stopcodon_usage_per_operon_position[position]
            tag_count = position_stopcodon_counts.get('TAG', 0)
            taa_count = position_stopcodon_counts.get('TAA', 0)
            tga_count = position_stopcodon_counts.get('TGA', 0)
            print position, taa_count, tag_count, tga_count
            write_handle.write('{0}\t{1}\t{2}\t{3}\n'.format(position, tag_count, taa_count, tga_count))

    print 'Done!'

if __name__ == '__main__':
    __main__()
