#!/usr/bin/env python
"""Module to run both concatemer_tree & compare_taxa."""

from divergence import parse_options
import os.path
import shutil
import sys
import tempfile

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


def main(args):
    """Main function called when run from command line or as part of pipeline."""

    #Input:
    #Taxon A
    #Taxon B
    #Retained orthologs

    #Output:
    #Orthologs per genome post-recombination filter
    #Concatemer post-recombination filter
    #Phylogenetic tree (now that branch lengths might have changed)

    #Tasks:
    #Run concatenate & deduce taxa
    #Run compare taxa with above produced post-filtering Taxon A & B

    usage = """
Usage: post_recombination_filter.py
--retained-orthologs-zip=FILE    archive of orthologous genes in FASTA format
--unfiltered-taxon-a=FILE        genome IDs for taxon A as deduced from phylogenetic tree of unfiltered concatemers
--unfiltered-taxon-b=FILE        genome IDs for taxon B as deduced from phylogenetic tree of unfiltered concatemers

--orthologs-per-genome=FILE      destination file path for orthologs split out per genome, based on the retained.zip
--concatemer=FILE                destination file path for super-concatemer of all genomes
--tree=FILE                      destination file path for tree visualization
"""
    options = ['retained-orthologs-zip', 'unfiltered-taxon-a', 'unfiltered-taxon-b',
               'orthologs-per-genome', 'concatemer', 'tree']
    retained_zip, unfiltered_a_file, unfiltered_b_file, \
    target_orth_per_genome, target_concat_file, target_tree = \
        parse_options(usage, options, args)

    #Run this in a temporary folder
    run_dir = tempfile.mkdtemp(prefix='post_recombination_filter_')

    #Create temporary files for the filtered taxon files
    filtered_a_file = os.path.join(run_dir, 'filtered_taxon_a.tsv')
    filtered_b_file = os.path.join(run_dir, 'filtered_taxon_b.tsv')

    #Run concatemer tree
    from concatemer_tree import main as ctree_main
    ctree_args = ['--orthologs-zip', retained_zip,
                  '--coding-regions', target_orth_per_genome,
                  '--concatemer', target_concat_file,
                  '--taxon-a', filtered_a_file,
                  '--taxon-b', filtered_b_file,
                  '--tree', target_tree]
    ctree_main(ctree_args)

    #Run compare taxa
    from compare_taxa import main as ctaxa_main
    ctaxa_args = ['--unfiltered-taxon-a', unfiltered_a_file,
                  '--unfiltered-taxon-b', unfiltered_b_file,
                  '--filtered-taxon-a', filtered_a_file,
                  '--filtered-taxon-b', filtered_b_file]
    ctaxa_main(ctaxa_args)

    #Remove unused files to free disk space
    shutil.rmtree(run_dir)

if __name__ == '__main__':
    main(sys.argv[1:])
