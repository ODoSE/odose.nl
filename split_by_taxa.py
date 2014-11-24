#!/usr/bin/env python
"""Module to split sequence files containing records from two taxa into separate sequence files per taxon."""

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from divergence import parse_options, extract_archive_of_files, create_directory, create_archive_of_files
import logging as log
import os.path
import re
import shutil
import sys
import tempfile

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


def split_alignment_by_taxa(run_dir, ortholog_files, (genome_ids_a, prefix_a), (genome_ids_b, prefix_b)):
    """Separate multiple sequence alignments with sequences from two taxa out into separate files per taxon."""
    # Collections to hold split files
    taxon_a_files = []
    taxon_b_files = []

    # For each alignment create separate alignments for clade A & clade B genomes
    for ortholog_file in ortholog_files:
        # Determine input file base name
        base_name, extension = os.path.splitext(os.path.split(ortholog_file)[1])

        # Separate alignment according to which taxon the genome_ids belong to
        alignment = AlignIO.read(ortholog_file, 'fasta')
        alignment_a = MultipleSeqAlignment(seqr for seqr in alignment if seqr.id.split('|')[0] in genome_ids_a)
        alignment_b = MultipleSeqAlignment(seqr for seqr in alignment if seqr.id.split('|')[0] in genome_ids_b)

        # Build up target output files
        taxon_a_file = os.path.join(run_dir, '{0}.{1}{2}'. format(base_name, prefix_a, extension))
        taxon_b_file = os.path.join(run_dir, '{0}.{1}{2}'. format(base_name, prefix_b, extension))

        # Actually write out sub alignments
        AlignIO.write(alignment_a, taxon_a_file, 'fasta')
        AlignIO.write(alignment_b, taxon_b_file, 'fasta')

        # Append the written files to the correct collections of files
        taxon_a_files.append(taxon_a_file)
        taxon_b_files.append(taxon_b_file)

    # Return collection of files
    return taxon_a_files, taxon_b_files


def _common_prefix(names, fallback=None):
    """Return common prefix from list of organism names with non alphanumeric characters removed, or fallback value
    if less than three characters are retained."""
    prefix = os.path.commonprefix(names).strip()
    clean_prefix = re.sub('[\W]', '', prefix)
    if len(clean_prefix) < 3:
        return fallback
    return clean_prefix


def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: split_by_taxa.py
--genomes-a=FILE        file with genome GenBank Project ID and Organism name on each line for taxon A
--genomes-b=FILE        file with genome GenBank Project ID and Organism name on each line for taxon B
--orthologs-zip=FILE    archive of aligned & trimmed single copy orthologous (SICO) genes
--taxon-a-zip=FILE      destination file path for archive of SICO genes belonging to taxon A
--taxon-b-zip=FILE      destination file path for archive of SICO genes belonging to taxon B
"""
    options = ['genomes-a', 'genomes-b', 'orthologs-zip', 'taxon-a-zip', 'taxon-b-zip']
    genome_a_ids_file, genome_b_ids_file, orthologs_zip, taxon_a_zip, taxon_b_zip = parse_options(usage, options, args)

    # Parse file containing RefSeq project IDs to extract RefSeq project IDs
    with open(genome_a_ids_file) as read_handle:
        lines = [line.split('\t') for line in read_handle]
        genome_ids_a = [line[0] for line in lines]
        common_prefix_a = _common_prefix([line[1] for line in lines], 'taxon_a')
    with open(genome_b_ids_file) as read_handle:
        lines = [line.split('\t') for line in read_handle]
        genome_ids_b = [line[0] for line in lines]
        common_prefix_b = _common_prefix([line[1] for line in lines], 'taxon_b')

    # Create run_dir to hold files related to this run
    run_dir = tempfile.mkdtemp(prefix='split_by_taxa_')

    # Extract files from zip archive
    ortholog_files = extract_archive_of_files(orthologs_zip, create_directory('alignments', inside_dir=run_dir))

    # Actually split alignments per taxon
    taxon_a_files, taxon_b_files = split_alignment_by_taxa(run_dir, ortholog_files,
                                                           (genome_ids_a, common_prefix_a),
                                                           (genome_ids_b, common_prefix_b))

    # Write the produced files to command line argument filenames
    create_archive_of_files(taxon_a_zip, taxon_a_files)
    create_archive_of_files(taxon_b_zip, taxon_b_files)

    # Remove unused files to free disk space
    shutil.rmtree(run_dir)

    # Exit after a comforting log message
    log.info("Produced: \n%s\n%s", taxon_a_zip, taxon_b_zip)
    return taxon_a_zip, taxon_b_zip

if __name__ == '__main__':
    main(sys.argv[1:])
