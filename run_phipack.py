#!/usr/bin/env python
"""Module to filter orthologs when recombination is found through PhiPack."""

from __future__ import division
from Bio import SeqIO
from shared import create_directory, extract_archive_of_files, parse_options, get_most_recent_gene_name, \
    find_cogs_in_sequence_records
from select_taxa import select_genomes_by_ids
from versions import PHIPACK
from subprocess import check_call, CalledProcessError
import logging as log
import os.path
import re
import shutil
import sys
import tempfile

__author__ = "Tim te Beek"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


def _phipack_for_all_orthologs(run_dir, aligned_files, stats_file):
    """Filter aligned fasta files where there is evidence of recombination when inspecting PhiPack values.
    Return two collections of aligned files, the first without recombination, the second with recombination."""

    log.info('Running PhiPack for %i orthologs to find recombination', len(aligned_files))

    # Create separate directory for phipack related values
    phipack_dir = create_directory('phipack', inside_dir=run_dir)

    with open(stats_file, mode='w') as write_handle:
        write_handle.write('\t'.join(['Ortholog',
                                      'Informative sites',
                                      'Phi',
                                      'Max Chi^2',
                                      'NSS',
                                      'COGs',
                                      'Product']) + '\n')

        # Retrieve unique genomes from first ortholog file
        genome_ids = set(fasta_record.id.split('|')[0] for fasta_record in SeqIO.parse(aligned_files[0], 'fasta'))
        genome_dicts = select_genomes_by_ids(genome_ids).values()

        # Assign ortholog files to the correct collection based on whether they show recombination
        for ortholog_file in aligned_files:
            orth_name = os.path.split(ortholog_file)[1].split('.')[0]

            # Parse tree file to ensure all genome_ids_a & genome_ids_b group together in the tree
            phipack_values = run_phipack(phipack_dir, ortholog_file)

            # Write PhiPack values to line
            write_handle.write('{0}\t{1[PhiPack sites]}\t{1[Phi]}\t{1[Max Chi^2]}\t{1[NSS]}'.format(orth_name,
                                                                                                    phipack_values))

            # Parse sequence records again, but now to retrieve cogs and products
            seq_records = list(SeqIO.parse(ortholog_file, 'fasta'))
            # COGs
            cogs = find_cogs_in_sequence_records(seq_records)
            write_handle.write('\t' + ','.join(cogs))
            # Product
            product = get_most_recent_gene_name(genome_dicts, seq_records)
            write_handle.write('\t' + product)

            # End line
            write_handle.write('\n')

    # Nothing to return, the stats_file is the product


def run_phipack(phipack_dir, dna_file):
    """Run PhiPack and return the number of informative sites, PHI, Max Chi^2 and NSS."""
    # Create directory for PhiPack to run in, so files get created there
    orth_name = os.path.split(dna_file)[1].split('.')[0]
    rundir = create_directory(orth_name, inside_dir=phipack_dir)

    # Build up list of commands
    command = PHIPACK, '-f', dna_file, '-o'  # Output NSS & Max Chi^2
    try:
        check_call(command, cwd=rundir, stdout=open('/dev/null', mode='w'))
    except CalledProcessError as err:
        log.warn('Error running PhiPack for %s:\n%s', dna_file, err)
        return {'PhiPack sites': None, 'Phi': None, 'Max Chi^2': None, 'NSS': None}

    # Retrieve output log file contents
    logfile = os.path.join(rundir, 'Phi.log')
    with open(logfile) as read_handle:
        contents = ''.join(read_handle)

    # Parse standard output to retrieve values for # sites, Phi, Chi^2 max & NSS
    # Found 103 informative sites.
    # PHI (Normal):        9.04e-01
    # Max Chi^2:           6.60e-01  (1000 permutations)
    # NSS:                 6.31e-01  (1000 permutations)
    sites = int(re.search('Found ([0-9]+) informative sites.', contents).group(1))
    raw_phi = re.search('PHI \(Normal\):\s+(.*)', contents).group(1)
    phi = float(raw_phi) if raw_phi != '--' else None
    chi = float(re.search('Max Chi\^2:\s+(.*)\s+\(1000 permutations\)', contents).group(1))
    nss = float(re.search('NSS:\s+(.*)\s+\(1000 permutations\)', contents).group(1))
    return {'PhiPack sites': sites, 'Phi': phi, 'Max Chi^2': chi, 'NSS': nss}


def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: run_phipack.py
--orthologs-zip=FILE     archive of orthologous genes in FASTA format
--stats-file=FILE        destination file path for values found through PhiPack for each ortholog
"""
    options = ('orthologs-zip', 'stats-file')
    orthologs_zip, stats_file = parse_options(usage, options, args)

    # Run filtering in a temporary folder, to prevent interference from simultaneous runs
    run_dir = tempfile.mkdtemp(prefix='run_phipack_')

    # Extract files from zip archive
    extraction_dir = create_directory('extracted_orthologs', inside_dir=run_dir)
    ortholog_files = extract_archive_of_files(orthologs_zip, extraction_dir)

    # Find recombination in all ortholog_files
    _phipack_for_all_orthologs(run_dir, ortholog_files, stats_file)

    # Remove unused files to free disk space
    shutil.rmtree(run_dir)

    # Exit after a comforting log message
    log.info('Produced:\n%s', stats_file)

if __name__ == '__main__':
    main(sys.argv[1:])
