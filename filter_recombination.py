#!/usr/bin/env python
"""Module to filter out orthologous genes with evidence of recombination between related taxa through inspecting 
phylogenetic trees."""

from Bio import AlignIO
from divergence import create_directory
from subprocess import Popen, PIPE, STDOUT
import os.path
import shutil
import tempfile

def filter_recombination_orthologs(gids_a, gids_b, aligned_files):
    """Filter aligned fasta files where there is evidence of recombination when inspecting phylogenetic trees."""
    run_dir = tempfile.mkdtemp(prefix = 'filter_recombination_')

    non_rec = [sico_file for sico_file in aligned_files if not _recombined_sico(run_dir, gids_a, gids_b, sico_file)]

    shutil.rmtree(run_dir)

    return non_rec


def _recombined_sico(run_dir, genome_ids_a, genome_ids_b, aligned_file):
    """Look for evidence of recombination between two taxa and return True if it is found, False if it is not found."""

    #Read alignment file
    alignment = AlignIO.read(aligned_file, 'fasta')

    #Create distance file
    distance_file = _run_dna_dist(run_dir, alignment)

    #Create tree file
    tree_file = _run_neighbor(run_dir, distance_file)

    #Parse tree file to ensure all genome_ids_a & genome_ids_b group together in the tree
    recombination_found = _find_recombination(genome_ids_a, genome_ids_b, tree_file)

    return recombination_found

def _run_dna_dist(run_dir, alignment):
    """Run dnadist to calculate distances between individual strains in a distance matrix, as input for neighbor."""
    #Run calculations inside a directory
    dnadist_dir = create_directory('dnadist/', inside_dir = run_dir)

    #Write input file for dnadist according to spec
    nr_of_species = len(alignment)
    nr_of_sites = len(alignment[0])
    infile = os.path.join(dnadist_dir, 'infile')
    with open(infile, mode = 'w') as write_handle:
        write_handle.write('   {0}   {2}'.format(nr_of_species, nr_of_sites))

        for seq_record in alignment:
            name = seq_record.id.split('|')[0]
            write_handle.write('{0:10}{1}'.format(name, seq_record.seq))

    #Actually run the dnadist program in the correct directory, and send input to it for the first prompt
    process = Popen('/path/to/dnadist', cwd = dnadist_dir, stdout = PIPE, stderr = STDOUT)
    process.communicate(input = 'Y\n')

    #Retrieve outputfile
    outfile = os.path.join(dnadist_dir, 'outfile')
    assert os.path.exists(outfile) and 0 < os.path.getsize(outfile), outfile + ' should exist with some content now'
    return outfile

def _run_neighbor(run_dir, distance_file):
    """Run neighbor to generate a tree of the distances in the distance file, and return the generated tree file."""
    neighbor_dir = create_directory('neighbor', inside_dir = run_dir)

    shutil.move(distance_file, os.path.join(neighbor_dir, 'infile'))

    #TODO Run neighbor

    tree_file = os.path.join(neighbor_dir, 'outtree')
    return tree_file


def _find_recombination(genome_ids_a, genome_ids_b, tree_file):
    """Look for evidence of recombination by seeing if all genomes of the separate taxa group together in the tree."""

    return False

if __name__ == '__main__':
    pass
