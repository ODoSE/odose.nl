#!/usr/bin/env python
"""Module to split multiple sequence alignments of two taxa into separate MSAs per taxon."""

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from divergence import parse_options, extract_archive_of_files, create_directory, create_archive_of_files
import logging as log
import os.path
import shutil
import sys
import tempfile

def split_alignment_by_taxon(run_dir, genome_ids_a, genome_ids_b, sico_files):
    """Separate mixed multiple sequence alignments with sequences from two taxa out into separate files per taxon."""
    #Collections to hold split files
    taxon_a_files = []
    taxon_b_files = []

    #Directories to hold split files
    taxon_a_dir = create_directory('taxon_a', inside_dir = run_dir)
    taxon_b_dir = create_directory('taxon_b', inside_dir = run_dir)

    #For each alignment create separate alignments for clade A & clade B genomes 
    for sico_file in sico_files:
        #Determine input file base name
        base_name = os.path.split(os.path.splitext(sico_file)[0])[1]

        #Separate alignment according to which taxon the genome_ids belong to
        alignment = AlignIO.read(sico_file, 'fasta')
        alignment_a = MultipleSeqAlignment(seqr for seqr in alignment if seqr.id.split('|')[0] in genome_ids_a)
        alignment_b = MultipleSeqAlignment(seqr for seqr in alignment if seqr.id.split('|')[0] in genome_ids_b)

        #Build up target output files
        taxon_a_file = os.path.join(taxon_a_dir, base_name + '_taxon_A.fasta')
        taxon_b_file = os.path.join(taxon_b_dir, base_name + '_taxon_B.fasta')

        #Actually write out sub alignments
        AlignIO.write(alignment_a, taxon_a_file, 'fasta')
        AlignIO.write(alignment_b, taxon_b_file, 'fasta')

        #Append the written files to the correct collections of files
        taxon_a_files.append(taxon_a_file)
        taxon_b_files.append(taxon_b_file)

    #Return collection of files
    return taxon_a_files, taxon_b_files

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: split_by_taxon.py 
--genomes-a=FILE     file with RefSeq id from complete genomes table on each line for taxon A
--genomes-b=FILE     file with RefSeq id from complete genomes table on each line for taxon B
--sico-zip=FILE      archive of aligned & trimmed single copy orthologous (SICO) genes
--taxon-a-zip=FILE    destination file path for archive of SICO genes belonging to taxon A
--taxon-b-zip=FILE    destination file path for archive of SICO genes belonging to taxon B
"""
    options = ['genomes-a', 'genomes-b', 'sico-zip', 'codeml-zip', 'dnds-stats']
    genome_a_ids_file, genome_b_ids_file, sico_zip, taxon_a_zip, taxon_b_zip = parse_options(usage, options, args)

    #Parse file containing RefSeq project IDs to extract RefSeq project IDs
    with open(genome_a_ids_file) as read_handle:
        genome_ids_a = [line.split()[0] for line in read_handle]
    with open(genome_b_ids_file) as read_handle:
        genome_ids_b = [line.split()[0] for line in read_handle]

    #Create run_dir to hold files related to this run
    run_dir = tempfile.mkdtemp(prefix = 'split_by_taxon_')

    #Extract files from zip archive
    sico_files = extract_archive_of_files(sico_zip, create_directory('alignments', inside_dir = run_dir))

    #Actually split alignments per taxon
    taxon_a_files, taxon_b_files = split_alignment_by_taxon(run_dir, genome_ids_a, genome_ids_b, sico_files)

    #Write the produced files to command line argument filenames
    create_archive_of_files(taxon_a_zip, taxon_a_files)
    create_archive_of_files(taxon_b_zip, taxon_b_files)

    #Remove unused files to free disk space 
    shutil.rmtree(run_dir)

    #Exit after a comforting log message
    log.info("Produced: \n%s\n%s", taxon_a_zip, taxon_b_zip)
    return taxon_a_zip, taxon_b_zip

if __name__ == '__main__':
    main(sys.argv[1:])
