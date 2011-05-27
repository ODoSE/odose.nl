#!/usr/bin/env python
"""Module to filter out orthologous genes with evidence of recombination between related taxa through inspecting 
phylogenetic trees."""

from Bio import AlignIO, Phylo
from divergence import create_directory
from subprocess import Popen, PIPE, STDOUT
import logging as log
import os.path
import shutil
import tempfile

def filter_recombined_orthologs(genome_ids_a, genome_ids_b, aligned_files):
    """Filter aligned fasta files where there is evidence of recombination when inspecting phylogenetic trees. 
    Return two collections of aligned files, the first without recombination, the second with recombination."""
    #Create temporary directory for run files
    run_dir = tempfile.mkdtemp(prefix = 'filter_recombination_')

    #Collections to hold both non recombination files & files showing recombination 
    non_recomb = []
    recombined = []

    #Assign ortholog files to the correct collection based on wether they show recombination 
    for ortholog_file in aligned_files:
        #Determine input file base name to create an ortholog run specific directory
        base_name = os.path.split(os.path.splitext(ortholog_file)[0])[1]
        ortholog_dir = create_directory(base_name, inside_dir = run_dir)

        #Create distance file
        distance_file = _run_dna_dist(ortholog_dir, ortholog_file)

        #Create tree file
        tree_file = _run_neighbor(ortholog_dir, distance_file)

        #Parse tree file to ensure all genome_ids_a & genome_ids_b group together in the tree
        if _find_recombination(genome_ids_a, genome_ids_b, tree_file):
            recombined.append(ortholog_file)
            log.info('%s showed evidence of recombination', base_name)
        else:
            non_recomb.append(ortholog_file)

    log.info('Recombination found in %i out of %i orthologs', len(recombined), len(aligned_files))

    #Remove temporary directory to clean up disk space
    shutil.rmtree(run_dir)

    return non_recomb, recombined

DNADIST = '/projects/divergence/software/phylip-3.69/exe/dnadist'

def _run_dna_dist(run_dir, aligned_file):
    """Run dnadist to calculate distances between individual strains in a distance matrix, as input for neighbor."""
    #Run calculations inside a directory
    dnadist_dir = create_directory('dnadist/', inside_dir = run_dir)

    #Read alignment file
    alignment = AlignIO.read(aligned_file, 'fasta')

    #Convert alignment in to proper input file for dnadist according to specification
    nr_of_species = len(alignment)
    nr_of_sites = len(alignment[0])
    infile = os.path.join(dnadist_dir, 'infile')
    with open(infile, mode = 'w') as write_handle:
        write_handle.write('   {0}   {1}\n'.format(nr_of_species, nr_of_sites))

        for seq_record in alignment:
            name = seq_record.id.split('|')[0]
            write_handle.write('{0:10}{1}\n'.format(name, seq_record.seq))

    #Actually run the dnadist program in the correct directory, and send input to it for the first prompt
    log.debug('Executing: %s in %s', DNADIST, dnadist_dir)
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

    #Copy outfile from dnadist to infile inside neighbor_dir
    shutil.copy(distance_file, os.path.join(neighbor_dir, 'infile'))

    #Actually run neighbor
    log.debug('Executing: %s in %s', NEIGHBOR, neighbor_dir)
    process = Popen(NEIGHBOR, cwd = neighbor_dir, stdin = PIPE, stdout = PIPE, stderr = STDOUT)
    process.communicate(input = 'N\nY\n')

    #Retrieve newick tree file
    treefile = os.path.join(neighbor_dir, 'outtree')
    assert os.path.exists(treefile) and 0 < os.path.getsize(treefile), treefile + ' should exist with some content now'
    return treefile

def _find_recombination(genome_ids_a, genome_ids_b, tree_file):
    """Look for evidence of recombination by seeing if all genomes of the separate taxa group together in the tree."""

    #Sample tree: (((((59245:0.00000,58803:0.00000):0.00000,59391:0.00000):0.01222,
    #((58531:0.00000,57915:0.00000):0.00182,(58623:0.00091,59379:0.00091):0.00091):0.01040):0.00229,
    #((59383:0.00457,58395:0.00457):0.00601,58783:0.01059):0.00393):0.06715,(58191:0.02967,(((59431:0.00367,
    #(58973:0.00000,58831:0.00000):0.00367):0.00117,(((58917:0.00000,59247:0.00000):0.00000,59249:0.00000):0.00367,
    #(59269:0.00000,58201:0.00000):0.00367):0.00117):0.00172,58017:0.00655):0.02311):0.05199)

    #Parse tree using BioPython, which wrongly interprets the RefSeq IDs as confidence scores, but that'll do for now.
    phylo_tree = Phylo.read(tree_file, 'newick')

    #Of the full tree retrieve the clades from the root clade, expecting exactly two distinct clades after UPGMA
    clades = phylo_tree.clade.clades
    assert len(clades) == 2

    #Get all the leafs for the above two clades in a similar format to the genome_ids
    clade_one = set(str(int(leaf.confidence)) for leaf in clades[0].get_terminals())
    clade_two = set(str(int(leaf.confidence)) for leaf in clades[1].get_terminals())

    #Use first genome of clade A to determine which collections should match with one another
    if genome_ids_a[0] in clade_one:
        #We'll declare to have found recombination when the taxa identified through the tree do not match the user taxa.
        recombination_found = set(genome_ids_a) != clade_one or set(genome_ids_b) != clade_two
    else:
        assert genome_ids_a[0] in clade_two
        recombination_found = set(genome_ids_a) != clade_one or set(genome_ids_b) != clade_two

    return recombination_found

if __name__ == '__main__':
    pass
