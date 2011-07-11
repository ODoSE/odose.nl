'''
Created on Jul 7, 2011

@author: tbeek

Hypothesis: There's a bias in the usage of stopcodons within operons. 
'''

from __future__ import division
from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_id
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

        if genome['Super Kingdom'] != 'Bacteria':
            continue

        #Skip some known troublesome cases
        if genome['RefSeq project ID'] in ('59133', '62947', '57853'):
            continue

        firstname = genome['Organism Name'].split()[0]

        #Skip genome if organism firstname was previously selected
        if firstname in previous_names:
            continue

        #Yield organism
        logging.info('Selected:\t%s', genome['Organism Name'])
        previous_names.add(firstname)
        yield genome

#Using the standard NCBI Bacterial, Archaeal and Plant Plastid Code translation table (11).
BACTERIAL_CODON_TABLE = unambiguous_dna_by_id.get(11)

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

def _get_stopcodon_usage(organism, genome_operons):
    """Get stopcodon usage per stopcodon in operon, and show difference between overall usage and as first/last."""
    #Some initial printed conclusions
    nr_operons = len(genome_operons)
    if nr_operons == 0:
        #Most probably no operons when the genome was a plasmid
        #print >> sys.stderr, 'No operons found'
        return
    nr_genes = sum(len(operon) for operon in genome_operons)
    if nr_genes == 0:
        print >> sys.stderr, 'No genes found'
        return

    #Extract only stopcodons
    stopcodons = [cds[2] for operon in genome_operons for cds in operon]

    #Get neutral usage statistics
    stopcodon_usage = dict((codon, stopcodons.count(codon)) for codon in set(stopcodons))
    assert sum(times for times in stopcodon_usage.values()) == nr_genes
    #Add percentages of total
    stopcodon_perc = dict((codon, times / nr_genes * 100) for codon, times in stopcodon_usage.iteritems())

    #Now determine the stopcodon usage for the first genes each operon
    first_stopcodons = [operon[0][2] for operon in genome_operons]
    first_stopcodon_perc = dict((codon2, times / nr_operons * 100) for codon2, times in
                                ((codon, first_stopcodons.count(codon)) for codon in set(first_stopcodons)))

    #Now determine the stopcodon usage for the first genes each operon
    last_stopcodons = [operon[-1][2] for operon in genome_operons]
    last_stopcodon_perc = dict((codon2, times / nr_operons * 100) for codon2, times in
                               ((codon, last_stopcodons.count(codon)) for codon in set(last_stopcodons)))

    for codon, perc in stopcodon_perc.iteritems():
        yield (organism,
               codon,
               nr_operons,
               nr_genes,
               nr_genes / nr_operons,
               stopcodon_usage[codon],
               stopcodon_perc[codon],
               first_stopcodon_perc.get(codon, 0) - perc,
               last_stopcodon_perc.get(codon, 0) - perc)

def __main__():
    """Main method called when run from terminal."""
    if True:
        #Select random genomes
        genomes = _sample_genomes(200)

        #Download genbank files for genomes in the background
        from multiprocessing import Pool
        pool = Pool()
        ft_genome_files = [pool.apply_async(download_genome_files, (genome,)) for genome in genomes]
        genome_files = (future.get() for future in ft_genome_files)
        genbank_files = (pairs[1] for files in genome_files for pairs in files)
    else:
        #Read current genbank from cache directory
        import glob
        genbank_files = list(glob.iglob('/data/dev/workspace-python/lib-divergence/divergence-cache/refseq/*/*.gbk'))

    #Extract coding sequences from genome files
    gbk_records = (SeqIO.read(genbank_file, 'genbank') for genbank_file in genbank_files)
    cds_per_genome = ((gbk_record.annotations['organism'], _extract_coding_sequences(gbk_record))
                                   for gbk_record in gbk_records)

    #Group coding sequences from each genome into operons
    operons_per_genome = ((organism, list(_group_cds_by_operons(genome_cds)))
                          for organism, genome_cds in cds_per_genome)

    #Per genome, calculate overall stopcodon usage
    with open('stats.tsv', mode = 'w', buffering = 1) as write_handle:
        header = '\t'.join(('Organism',
                            'Codon',
                            'Operons',
                            'Genes',
                            'Genes per Operon',
                            'Times used overall',
                            'Fraction of stopcodons',
                            'Percentage +/- as first',
                            'Percentage +/- as last'))
        print header
        write_handle.write(header + ' \n')

        taa_tuples = []
        tag_tuples = []
        tga_tuples = []

        usage_tuples = (_get_stopcodon_usage(organism, genm_operons) for organism, genm_operons in operons_per_genome)
        codongetter = itemgetter(1)
        for usage_tuple in chain.from_iterable(usage_tuples):
            if codongetter(usage_tuple) == 'TAA':
                taa_tuples.append(usage_tuple)
            elif codongetter(usage_tuple) == 'TAG':
                tag_tuples.append(usage_tuple)
            elif codongetter(usage_tuple) == 'TGA':
                tga_tuples.append(usage_tuple)
            line = '\t'.join(str(item) for item in usage_tuple)
            print line
            write_handle.write(line + '\n')

        #Write out averages
        write_handle.write(header + ' \n')
        taa_averages = (sum(usage_tuple[idx] for usage_tuple in taa_tuples) / len(taa_tuples) for idx in range(2, 9))
        write_handle.write('Average\tTAA\t' + '\t'.join(str(average) for average in taa_averages) + ' \n')
        tag_averages = (sum(usage_tuple[idx] for usage_tuple in tag_tuples) / len(tag_tuples) for idx in range(2, 9))
        write_handle.write('Average\tTAG\t' + '\t'.join(str(average) for average in tag_averages) + ' \n')
        tga_averages = (sum(usage_tuple[idx] for usage_tuple in tga_tuples) / len(tga_tuples) for idx in range(2, 9))
        write_handle.write('Average\tTGA\t' + '\t'.join(str(average) for average in tga_averages) + ' \n')

    print 'Done!'

if __name__ == '__main__':
    __main__()
