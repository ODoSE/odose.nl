#!/usr/bin/env python
"""Module to extract from DNA files orthologous sequences according to groups found through OrthoMCL step."""

from __future__ import division
from Bio import SeqIO
from divergence import create_directory, extract_archive_of_files, create_archive_of_files, parse_options
from itertools import chain
import logging as log
import os
import shutil
import sys
import tempfile

def _produce_heatmap(genomes, shared_single_copy, shared_multi_copy, accessory_genes):

    print shared_single_copy

    sico_counts = ((str(len(genome_proteins.get(genome, []))) for genome in genomes) for genome_proteins in shared_single_copy)
    muco_counts = (tuple(str(len(genome_proteins.get(genome, []))) for genome in genomes) for genome_proteins in shared_multi_copy)
    accessory_counts = (tuple(str(len(genome_proteins.get(genome, []))) for genome in genomes) for genome_proteins in accessory_genes)

    print '\t'.join(genomes)
    for some_count in sico_counts:
        print '\t'.join(some_count)
    for some_count in sorted(muco_counts):
        print '\t'.join(some_count)
    for some_count in sorted(accessory_counts):
        print '\t'.join(some_count)
    pass

def extract_orthologs(run_dir, genomes, dna_files, groups_file):
    """Extract DNA sequences for SICO, MUCO & partially shared orthologs to a single file per ortholog."""
    #Subdivide orthologs into groups
    shared_single_copy, shared_multi_copy, accessory_genes = _extract_shared_orthologs(genomes, groups_file)

    #Produce heatmap
    _produce_heatmap(genomes, shared_single_copy, shared_multi_copy, accessory_genes)

    #Extract fasta files per orthologs
    sico_files, muco_files, subset_files, nr_of_seqs = \
        _dna_file_per_sico(run_dir, dna_files, shared_single_copy, shared_multi_copy, accessory_genes)

    #Assertions
    if shared_single_copy:
        assert sico_files
    if shared_multi_copy:
        assert muco_files
    if accessory_genes:
        assert subset_files

    #Write statistics file
    stats_file = _write_statistics_file(run_dir, genomes, shared_single_copy, shared_multi_copy, accessory_genes, nr_of_seqs)

    return sico_files, muco_files, subset_files, stats_file

def _create_ortholog_dictionaries(groups_file):
    """Convert groups file into a list of ortholog dictionaries, which map project_id to their associated proteins."""
    #Sample line: 58017|YP_219088.1 58191|YP_001572431.1 59431|YP_002149136.1
    ortholog_proteins_per_genome = []
    with open(groups_file) as read_handle:
        for line in read_handle:
            remainder = line.split()
            proteins_per_genome = {}
            for ortholog in remainder:
                #Sample: 58017|YP_219088.1
                project_id, protein_id = ortholog.split('|')
                #Use dict().get(key, fallback_value) here to retrieve and assign valid array values for missing keys 
                proteins_per_genome[project_id] = proteins_per_genome.get(project_id, [])
                proteins_per_genome[project_id].append(protein_id)
            #Assign proteins per genome dictionary to orthologs per group as 
            ortholog_proteins_per_genome.append(proteins_per_genome)
    return ortholog_proteins_per_genome

def _extract_shared_orthologs(selected_genome_ids, groups_file):
    """Filter orthologs to retain shared single and multiple copy orthologs from the collection of genomes."""
    log.info('Extracting shared orthologs for %d genomes from %s', len(selected_genome_ids), groups_file)
    ortholog_proteins_per_genome = _create_ortholog_dictionaries(groups_file)

    #Group orthologs into the following categories
    shared_single_copy = []
    shared_multi_copy = []
    non_shared_orthologs = []
    for prot_per_genomes in ortholog_proteins_per_genome:
        #Assume both shared and single copy
        is_shared_genome = True
        is_single_copy = True

        #Test validity of above boolean statements against proteins per genome 
        for selected_genome in selected_genome_ids:
            if selected_genome not in prot_per_genomes:
                #If any selected genome does not have proteins in this ortholog, it is not shared (no need to go on)
                is_shared_genome = False
                break
            elif 1 < len(prot_per_genomes[selected_genome]):
                #If any selected genome has more than one protein in this ortholog, it is not a single copy ortholog
                is_single_copy = False

        #Based on the now validated above boolean statements, optionally add proteins per genome to shared collections
        if is_shared_genome:
            if is_single_copy:
                shared_single_copy.append(prot_per_genomes)
            else:
                shared_multi_copy.append(prot_per_genomes)
        else:
            non_shared_orthologs.append(prot_per_genomes)

    #Assert all orthologs are binned
    sum_length = len(shared_single_copy) + len(shared_multi_copy) + len(non_shared_orthologs)
    assert len(ortholog_proteins_per_genome) == sum_length, 'All orthologs should fall into any one of three groups'

    #Return three collections of dictionaries mapping project_id to proteins for all orthologs
    return shared_single_copy, shared_multi_copy, non_shared_orthologs

def _dna_file_per_sico(run_dir, dna_files, shared_single_copy, shared_multi_copy, non_shared):
    """Create fasta files with all sequences per ortholog."""
    #Delete & create directory to remove any previously existing SICO files
    sico_dir = create_directory('sico', inside_dir = run_dir)
    muco_dir = create_directory('muco', inside_dir = run_dir)
    subset_dir = create_directory('subset', inside_dir = run_dir)

    #Loop over DNA files to extract SICO genes from each genome to file per SICO
    sico_files = set()
    muco_files = set()
    subset_files = set()
    number_of_sequences = 0
    for dna_file in dna_files:
        log.info('Extracting orthologous genes from %s', dna_file)
        for record in SeqIO.parse(dna_file, 'fasta'):
            number_of_sequences += 1

            #Find record in each list of dictionaries, to append it to the corresponding ortholog files 
            sico_files.update(_write_record_to_ortholog_file(sico_dir, shared_single_copy, record))
            muco_files.update(_write_record_to_ortholog_file(muco_dir, shared_multi_copy, record))
            subset_files.update(_write_record_to_ortholog_file(subset_dir, non_shared, record))

    return sorted(sico_files), sorted(muco_files), sorted(subset_files), number_of_sequences

def _write_record_to_ortholog_file(directory, ortholog_dictionaries, record):
    """Find sequence record in list of ortholog dictionaries, to write the record to corresponding ortholog file."""
    #Sample header line:    >58191|NC_010067.1|YP_001569097.1|COG4948MR|core
    #Corresponding ortholog: {'58191': ['YP_001569097.1'], ...}
    split_header = record.id.split('|')
    project_id = split_header[0]
    protein_id = split_header[2]

    affected_ortholog_files = set()
    #Find the current sequence record in dictionaries of orthologs to append it to the right target file(s?)
    for number, ortholog in enumerate(ortholog_dictionaries):
        #Skip this ortholog if it does not contain project_id
        if project_id not in ortholog:
            continue

        #Write header & sequence to ortholog file if protein ID is included in protein IDs mapped to by project_id
        if protein_id in ortholog[project_id]:
            #Create filename for the current ortholog, to which previous or further sequence records might also be added
            sico_file = os.path.join(directory, 'ortholog_{0:06}.ffn'.format(number))
            affected_ortholog_files.add(sico_file)

            #Append to the MUCO files here to group the orthologs from various genomes in the same file
            with open(sico_file, mode = 'a') as write_handle:
                SeqIO.write(record, write_handle, 'fasta')

            #Can't rule out that a single protein is part of multiple orthologs, without intricate knowledge of OrthoMCL
            #Therefore continue looking for further occurrences of this protein in other orthologs
            continue

    #Return affected files, so a complete set of ortholog files can be build by caller
    return affected_ortholog_files

def _write_statistics_file(run_dir, genomes, shared_single_copy, shared_multi_copy, partially_shared, nr_of_seqs):
    """Write out file with some basic statistics about the genomes, orthologs and size of shared core genome."""
    #Some easy statistics about genomes and orthologs
    nr_shared_sico = len(shared_single_copy)
    nr_shared_muco = len(shared_multi_copy)
    nr_part_shared = len(partially_shared)
    nr_orthologs = nr_shared_sico + nr_shared_muco + nr_part_shared

    #Count the number of partially shared orthologs per number of genomes shared among
    orthologs_per_nr_of_genomes = {}
    for ortholog in partially_shared:
        nr_of_genomes = len(ortholog)
        orthologs_per_nr_of_genomes[nr_of_genomes] = orthologs_per_nr_of_genomes.get(nr_of_genomes, 0) + 1

    #Determine number of ORFans by deducting unique proteins identified as orthologs from total number of genes
    proteins = set(chain.from_iterable(prot for per_genome in shared_single_copy for prot in per_genome.values()))
    proteins.update(chain.from_iterable(prot for per_genome in shared_multi_copy for prot in per_genome.values()))
    proteins.update(chain.from_iterable(prot for per_genome in partially_shared for prot in per_genome.values()))
    nr_orfans = nr_of_seqs - len(proteins)

    stats_file = os.path.join(run_dir, 'extract-stats.txt')
    with open(stats_file, mode = 'w') as writer:
        #Write Genome & gene count statistics to file        
        writer.write('#{0:7}\tGenomes\n'.format(len(genomes)))
        writer.write('#{0:7}\tGenes\n'.format(nr_of_seqs))
        writer.write('#{0:7}\tORFan genes (no orthologs)\n\n'.format(nr_orfans))

        #Write Ortholog count statistics to file
        def perc(number):
            """Calculate a number as percentage of the number of the number of orthologs"""
            return number / nr_orthologs
        writer.write('# Distribution of orthologs over identified groups:\n')
        writer.write('{0:8}\t{1:8.2%}\tSingle-copy orthologs shared across all genomes\n'
                     .format(nr_shared_sico, perc(nr_shared_sico)))
        writer.write('{0:8}\t{1:8.2%}\tMultiple-copy orthologs shared across all genomes\n'
                     .format(nr_shared_muco, perc(nr_shared_muco)))
        #Print the number of orthologs found per nr of genomes
        for nr_of_genomes, nr_of_orthologs in orthologs_per_nr_of_genomes.iteritems():
            writer.write('{0:8}\t{1:8.2%}\t{2}-genome orthologs with single or multiple copies.\n'
                     .format(nr_of_orthologs, perc(nr_of_orthologs), nr_of_genomes))
        writer.write('{0:8}\t{1:8.2%}\tTotal number of orthologs\n'.format(nr_orthologs, perc(nr_orthologs)))

    assert os.path.isfile(stats_file) and 0 < os.path.getsize(stats_file), stats_file + ' should exist with content.'
    return stats_file

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: extract_orthologs.py 
--genomes=FILE       file with GenBank Project IDs from complete genomes table on each line 
--dna-zip=FILE       zip archive of extracted DNA files
--groups=FILE        file listing groups of orthologous proteins

--sico-zip=FILE      destination file path for archive of shared single copy orthologous (SICO) genes
--muco-zip=FILE      destination file path for archive of shared multiple copy orthologous genes
--subset-zip=FILE    destination file path for archive of variable copy orthologous genes shared for a subset only
--stats=FILE         destination file path for ortholog statistics file
"""
    options = ['genomes', 'dna-zip', 'groups', 'sico-zip', 'muco-zip', 'subset-zip', 'stats']
    genome_ids_file, dna_zip, groups_file, target_sico, target_muco, target_subset, target_stats_path = \
    parse_options(usage, options, args)

    #Parse file extract GenBank Project IDs
    with open(genome_ids_file) as read_handle:
        genomes = [line.split()[0] for line in read_handle if not line.startswith('#')]

    #Create temporary directory within which to extract orthologs
    run_dir = tempfile.mkdtemp(prefix = 'extract_orthologs_run_')

    #Extract files from zip archive
    temp_dir = create_directory('dna_files', inside_dir = run_dir)
    dna_files = extract_archive_of_files(dna_zip, temp_dir)

    #Actually run ortholog extraction
    sico_files, muco_files, subset_files, stats_file = extract_orthologs(run_dir, genomes, dna_files, groups_file)

    #Move produced files to command line specified output paths
    create_archive_of_files(target_sico, sico_files)
    create_archive_of_files(target_muco, muco_files)
    create_archive_of_files(target_subset, subset_files)
    shutil.move(stats_file, target_stats_path)

    #Remove unused files to free disk space 
    shutil.rmtree(run_dir)

    #Exit after a comforting log message
    log.info("Produced: \n%s\n%s\n%s\n%s", target_sico, target_muco, target_subset, target_stats_path)

if __name__ == '__main__':
    main(sys.argv[1:])
