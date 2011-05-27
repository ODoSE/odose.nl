#!/usr/bin/env python
"""Module to filter out orthologous genes with evidence of recombination between related taxa through inspecting 
phylogenetic trees."""

from Bio import AlignIO
from divergence import create_directory
from subprocess import Popen, PIPE, STDOUT
import os.path
import shutil
import tempfile
import logging as log

def filter_recombined_orthologs(gids_a, gids_b, aligned_files):
    """Filter aligned fasta files where there is evidence of recombination when inspecting phylogenetic trees."""
    run_dir = tempfile.mkdtemp(prefix = 'filter_recombination_')

    non_rec = [sico_file for sico_file in aligned_files if not _is_recombined(run_dir, gids_a, gids_b, sico_file)]

    shutil.rmtree(run_dir)

    return non_rec


def _is_recombined(run_dir, genome_ids_a, genome_ids_b, aligned_file):
    """Look for evidence of recombination between two taxa and return True if it is found, False if it is not found."""

    #Determine input file base name to create a run specific directory
    base_name = os.path.split(os.path.splitext(aligned_file)[0])[1]
    ortholog_dir = create_directory(base_name, inside_dir = run_dir)

    #Read alignment file
    alignment = AlignIO.read(aligned_file, 'fasta')

    #Create distance file
    distance_file = _run_dna_dist(ortholog_dir, alignment)

    #Create tree file
    tree_file = _run_neighbor(ortholog_dir, distance_file)

    #Parse tree file to ensure all genome_ids_a & genome_ids_b group together in the tree
    recombination_found = _find_recombination(genome_ids_a, genome_ids_b, tree_file)

    return recombination_found

DNADIST = '/projects/divergence/software/phylip-3.69/exe/dnadist'

def _run_dna_dist(run_dir, alignment):
    """Run dnadist to calculate distances between individual strains in a distance matrix, as input for neighbor."""
    #Run calculations inside a directory
    dnadist_dir = create_directory('dnadist/', inside_dir = run_dir)

    #Write input file for dnadist according to spec
    nr_of_species = len(alignment)
    nr_of_sites = len(alignment[0])
    infile = os.path.join(dnadist_dir, 'infile')
    with open(infile, mode = 'w') as write_handle:
        write_handle.write('   {0}   {1}\n'.format(nr_of_species, nr_of_sites))

        for seq_record in alignment:
            name = seq_record.id.split('|')[0]
            write_handle.write('{0:10}{1}\n'.format(name, seq_record.seq))

    #Actually run the dnadist program in the correct directory, and send input to it for the first prompt
    log.info('Executing: %s in %s', DNADIST, dnadist_dir)
    process = Popen(DNADIST, cwd = dnadist_dir, stdin = PIPE, stdout = PIPE, stderr = STDOUT)
    process.communicate(input = 'Y\n')

    #Retrieve outputfile
    outfile = os.path.join(dnadist_dir, 'outfile')
    assert os.path.exists(outfile) and 0 < os.path.getsize(outfile), outfile + ' should exist with some content now'
    return outfile

NEIGHBOR = '/projects/divergence/software/phylip-3.69/exe/neighbor'

def _run_neighbor(run_dir, distance_file):
    """Run neighbor to generate a tree of the distances in the distance file, and return the generated tree file."""
    neighbor_dir = create_directory('neighbor', inside_dir = run_dir)

    shutil.move(distance_file, os.path.join(neighbor_dir, 'infile'))

    #Actually run neighbor
    log.info('Executing: %s in %s', NEIGHBOR, neighbor_dir)
    process = Popen(NEIGHBOR, cwd = neighbor_dir, stdin = PIPE, stdout = PIPE, stderr = STDOUT)
    process.communicate(input = 'Y\n')

    #Retrieve newick tree file
    treefile = os.path.join(neighbor_dir, 'outtree')
    assert os.path.exists(treefile) and 0 < os.path.getsize(treefile), treefile + ' should exist with some content now'
    return treefile


def _find_recombination(genome_ids_a, genome_ids_b, tree_file):
    """Look for evidence of recombination by seeing if all genomes of the separate taxa group together in the tree."""

    with open(tree_file) as read_handle:
        tree = ''.join(line.strip() for line in read_handle)
        print tree


    return False

if __name__ == '__main__':
    pass
