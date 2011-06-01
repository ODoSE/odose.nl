#!/usr/bin/env python
"""Module to deduce taxa from genome concatemer using phylogenetic tree and the clades identified therein."""

from Bio import Phylo
from divergence import parse_options
from divergence.filter_orthologs import _run_dna_dist, _run_neighbor, _read_taxa_from_tree
from divergence.select_taxa import select_genomes_by_ids
import logging as log
import sys
import tempfile
import shutil

def visualize_tree(super_tree_file, id_to_name_map, ascii_tree):
    """Visualize the phylogenetic tree encoded in the Newick format super_tree_file, and write graphic to ascii_tree."""
    #Draw phylogenetic tree
    tree = Phylo.read(super_tree_file, 'newick')

    #BioPython misinterprets numerical leaf names as confidence scores: Fix this here
    for leaf in tree.get_terminals():
        refseq_id = str(int(leaf.confidence))
        leaf.name = '{0}\t{1}'.format(refseq_id, id_to_name_map[refseq_id])

    #Print ascii tree, as we can't get visualization to work properly using draw_graphviz
    with open(ascii_tree, mode = 'w') as write_handle:
        Phylo.draw_ascii(tree, file = write_handle, column_width = 120)

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: deduce_taxa_from_tree.py
--concatemer=FILE    super-concatemer of all genomes in FASTA format
--taxon-a=FILE       destination file path for genome IDs for taxon A
--taxon-b=FILE       destination file path for genome IDs for taxon B
--tree=FILE          destination file path for tree visualization
"""
    options = ['concatemer', 'taxon-a', 'taxon-b', 'tree']
    concatemer, target_taxon_a, target_taxon_b, target_tree = parse_options(usage, options, args)

    #Director to contain run specific files
    run_dir = tempfile.mkdtemp(prefix = 'deduce_taxa_')

    #Determine the taxa present in the super concatemer tree by building a phylogenetic tree from genome concatemer and
    #reading genome ids in the two largest clades.
    super_distance_file = _run_dna_dist(run_dir, concatemer)
    super_tree_file = _run_neighbor(run_dir, super_distance_file)
    genome_ids_a, genome_ids_b = _read_taxa_from_tree(super_tree_file)

    #Map RefSeq project IDs to Organism names
    id_to_name_map = dict((genome['RefSeq project ID'], genome['Organism Name']) for genome in \
                          select_genomes_by_ids(genome_ids_a + genome_ids_b))

    #Write RefSeq project ID and Organism Name to files
    with open(target_taxon_a, mode = 'w') as write_handle:
        for genome_id in genome_ids_a:
            write_handle.write('{id}\t{name}\n'.format(id = genome_id, name = id_to_name_map[genome_id]))
    with open(target_taxon_b, mode = 'w') as write_handle:
        for genome_id in genome_ids_b:
            write_handle.write('{id}\t{name}\n'.format(id = genome_id, name = id_to_name_map[genome_id]))

    #Visualize tree
    visualize_tree(super_tree_file, id_to_name_map, target_tree)

    #Clean up
    shutil.rmtree(run_dir)

    #Exit after a comforting log message
    log.info('Produced: \n%s\n%s\n%s', target_taxon_a, target_taxon_b, target_tree)
    return target_taxon_a, target_taxon_b, target_tree

if __name__ == '__main__':
    main(sys.argv[1:])
