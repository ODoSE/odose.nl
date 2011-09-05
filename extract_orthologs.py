#!/usr/bin/env python
###
# Part of the Adaptive Divergence through Direction of Selection workflow.
# Copyright (C) 2011  Tim te Beek <tim.te.beek@nbic.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###
"""Module to extract from DNA files orthologous sequences according to groups found through OrthoMCL step."""

from __future__ import division
from Bio import SeqIO
from divergence import create_directory, extract_archive_of_files, create_archive_of_files, parse_options
from divergence.filter_orthologs import find_cogs_in_sequence_records
from itertools import chain
import logging as log
import os
import shutil
import sys
import tempfile

def _produce_heatmap(genome_ids, sico_files, muco_files, accessory_files):
    """Produce heatmap of orthologs, and how many times ortholog ooccurs in genome, with the COGs added as well. """
    def _occurences_and_cogs(genome_ids, ortholog_files):
        """Generator that returns how many sequences exist per genome in each ortholog in order and which COGs occur."""
        for fasta_file in ortholog_files:
            records = tuple(SeqIO.parse(fasta_file, 'fasta'))
            ids = [record.id.split('|')[0] for record in records]
            count_per_id = [ids.count(genome_id) for genome_id in genome_ids]
            cogs = sorted(find_cogs_in_sequence_records(records))
            yield count_per_id, cogs

    heatmap = tempfile.mkstemp(suffix = '.tsv', prefix = 'genome_ortholog_heatmap_')[1]
    with open(heatmap, mode = 'w') as write_handle:
        #Write file header
        write_handle.write('\t'.join(genome_ids))
        write_handle.write('\tCOGs\n')

        #Write out sico
        for counts_per_id, cogs in sorted(_occurences_and_cogs(genome_ids, sico_files)):
            write_handle.write('\t'.join(str(occurrences) for occurrences in counts_per_id))
            write_handle.write('\t' + ','.join(cogs) + '\n')

        #Write out muco
        for counts_per_id, cogs in sorted(_occurences_and_cogs(genome_ids, muco_files)):
            write_handle.write('\t'.join(str(occurrences) for occurrences in counts_per_id))
            write_handle.write('\t' + ','.join(cogs) + '\n')

        #Write out accessory
        for counts_per_id, cogs in sorted(_occurences_and_cogs(genome_ids, accessory_files)):
            write_handle.write('\t'.join(str(occurrences) for occurrences in counts_per_id))
            write_handle.write('\t' + ','.join(cogs) + '\n')
    return heatmap

def extract_orthologs(run_dir, genomes, dna_files, groups_file, require_limiter = False):
    """Extract DNA sequences for SICO, MUCO & partially shared orthologs to a single file per ortholog."""
    #Subdivide orthologs into groups
    shared_single_copy, shared_multi_copy, accessory = _extract_shared_orthologs(genomes, groups_file, require_limiter)

    #Extract fasta files per orthologs
    sico_files, muco_files, accessory_files, nr_of_seqs = \
        _dna_file_per_sico(run_dir, dna_files, shared_single_copy, shared_multi_copy, accessory)

    #Produce heatmap
    heatmap_file = _produce_heatmap(genomes, sico_files, muco_files, accessory_files)

    #Assertions
    if shared_single_copy:
        assert sico_files
    if shared_multi_copy:
        assert muco_files
    if accessory:
        assert accessory_files

    #Write statistics file
    stats_file = _write_statistics_file(run_dir, genomes, shared_single_copy, shared_multi_copy, accessory, nr_of_seqs)

    return sico_files, muco_files, accessory_files, stats_file, heatmap_file

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

def _extract_shared_orthologs(selected_genome_ids, groups_file, require_limiter_presence = False):
    """Filter orthologs to retain shared single and multiple copy orthologs from the collection of genomes."""
    log.info('Extracting shared orthologs for %d genomes from %s', len(selected_genome_ids), groups_file)
    ortholog_proteins_per_genome = _create_ortholog_dictionaries(groups_file)

    #Group orthologs into the following categories
    shared_single_copy = []
    shared_multi_copy = []
    non_shared_orthologs = []
    for prot_per_genomes in ortholog_proteins_per_genome:
        #The ortholog is shared if all selected ids are present in the keys from the prot_per_genomes dictionary
        #While this still allows the selected_genome_ids to be a subset of all genome ids present in these orthologs
        is_shared_ortholog = all(selected_id in prot_per_genomes
                               for selected_id in selected_genome_ids)

        #If limiter presence is required, only mark orthologs as shared if they contain also limiter
        if require_limiter_presence:
            is_shared_ortholog = is_shared_ortholog and 'limiter' in prot_per_genomes

        #The ortholog is single copy if for all genomes in selected genome ids the number of proteins is exactly one
        #While this still allows for multiple copies in proteins not included in selected_genome_ids   
        is_single_copy = all(len(proteins) == 1
                             for genome_id, proteins in prot_per_genomes.iteritems()
                             if genome_id in selected_genome_ids)

        #Based on the now validated above boolean statements, optionally add proteins per genome to shared collections
        if is_shared_ortholog:
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
    nr_shared_sico_orth = len(shared_single_copy)
    nr_non_sico_orth = len(shared_multi_copy) + len(partially_shared)

    #Determine number of ORFans by deducting unique proteins identified as orthologs from total number of genes
    proteins = set(chain.from_iterable(prot for per_genome in shared_single_copy for prot in per_genome.values()))
    nr_sico_genes = len(proteins)
    proteins.update(chain.from_iterable(prot for per_genome in shared_multi_copy for prot in per_genome.values()))
    proteins.update(chain.from_iterable(prot for per_genome in partially_shared for prot in per_genome.values()))
    nr_non_sico_genes = len(proteins) - nr_sico_genes
    nr_orfans = nr_of_seqs - len(proteins)

    stats_file = os.path.join(run_dir, 'extract-stats.txt')
    with open(stats_file, mode = 'w') as writer:
        #Write Genome & gene count statistics to file        
        writer.write('#{0:6}\tGenomes\n'.format(len(genomes)))
        writer.write('{0:7}\tGenes\n'.format(nr_of_seqs))
        writer.write('{0:7}\tORFan genes (no orthologs)\n'.format(nr_orfans))
        writer.write('{0:7}\tShared single-copy orthologous genes in {1} orthologs\n'.format(nr_sico_genes,
                                                                                               nr_shared_sico_orth))
        writer.write('{0:7}\tOtherwise orthologous genes in {1} orthologs\n'.format(nr_non_sico_genes,
                                                                                      nr_non_sico_orth))

    assert os.path.isfile(stats_file) and 0 < os.path.getsize(stats_file), stats_file + ' should exist with content.'
    return stats_file

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: extract_orthologs.py 
--genomes=FILE       file with GenBank Project IDs from complete genomes table on each line 
--dna-zip=FILE       zip archive of extracted DNA files
--groups=FILE        file listing groups of orthologous proteins
--require-limiter    flag whether extracted core set of genomes should contain the limiter added in OrthoMCL [OPTIONAL] 

--sico-zip=FILE      destination file path for archive of shared single copy orthologous (SICO) genes
--muco-zip=FILE      destination file path for archive of shared multiple copy orthologous genes
--subset-zip=FILE    destination file path for archive of variable copy orthologous genes shared for a subset only
--stats=FILE         destination file path for ortholog statistics file
--heatmap=FILE       destination file path heatmap of orthologs and occurrences of ortholog per genome
"""
    options = ['genomes', 'dna-zip', 'groups', 'require-limiter?',
               'sico-zip', 'muco-zip', 'subset-zip', 'stats', 'heatmap']
    genome_ids_file, dna_zip, groups_file, require_limiter, \
    target_sico, target_muco, target_subset, target_stats_path, target_heat = \
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
    sico_files, muco_files, subset_files, stats_file, heatmap_file = extract_orthologs(run_dir, genomes, dna_files,
                                                                                       groups_file, require_limiter)

    #Move produced files to command line specified output paths
    create_archive_of_files(target_sico, sico_files)
    create_archive_of_files(target_muco, muco_files)
    create_archive_of_files(target_subset, subset_files)
    shutil.move(stats_file, target_stats_path)
    shutil.move(heatmap_file, target_heat)

    #Remove unused files to free disk space 
    shutil.rmtree(run_dir)

    #Exit after a comforting log message
    log.info("Produced: \n%s\n%s\n%s\n%s\n%s", target_sico, target_muco, target_subset, target_stats_path, target_heat)

if __name__ == '__main__':
    main(sys.argv[1:])
