#!/usr/bin/env python
"""Module for the select taxa step."""

import logging
from operator import itemgetter
import sys

from download_taxa_ncbi import download_genome_files
from load_prokaryotes import _parse_genomes_table
from shared import parse_options


__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


def select_genomes_by_ids(genome_ids):
    """Return list of genomes from complete genomes table whose Assembly Accession is in genome_ids."""
    # Loop over genomes and return any genomes whose GenBank Project ID is in genome_ids
    refseq_genomes = dict((genome['Assembly Accession'], genome) for genome in _parse_genomes_table())

    # Match genomes_ids to genomes
    matches = dict((queryid, refseq_genomes[queryid]) for queryid in genome_ids if queryid in refseq_genomes)

    # See if we matched all genomes, else log a warning
    for queryid in genome_ids:
        if queryid not in matches:
            logging.warning('Could not find genome with BioProject ID %s in complete genomes table', queryid)

    return matches


def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: select_taxa.py
--genomes=ID,...           optional comma-separated list of selected GenBank Project IDs from complete genomes table
--previous-file=FILE       optional previously or externally created GenBank Project IDs file whose genomes should be reselected
--require-protein-table    require protein table files to be present for all downloaded genomes
--genomes-file=FILE        destination path for file with selected genome IDs followed by Organism Name on each line
"""
    options = ['genomes=?', 'previous-file=?', 'require-protein-table?', 'genomes-file']
    genomes_line, previous_file, require_ptt, genomes_file = parse_options(usage, options, args)

    # Genome IDs selected by the user that refer to GenBank or RefSeq entries
    genome_ids = []

    # Split ids on comma
    if genomes_line:
        genome_ids.extend(val for val in genomes_line.split(',') if val)

    # Allow for input of previous or externally created genomes-file to rerun an analysis
    if previous_file:
        # Read previous GenBank Project IDs from previous_file, each on their own line
        with open(previous_file) as read_handle:
            genome_ids.extend(line.split()[0] for line in read_handle
                              # But skip external genomes as their IDs will fail to download
                              if 'external genome' not in line)

    # Assert each clade contains enough IDs
    maximum = 100
    # TODO Move this test to translate, where we can see how many translations succeeded + how many externals there are
    if maximum < len(genome_ids):
        logging.error('Expected between two and {0} selected genomes, but was {1}'.format(maximum, len(genome_ids)))
        sys.exit(1)

    # Retrieve genome dictionaries to get to Organism Name
    genomes = select_genomes_by_ids(genome_ids).values()
    genomes = sorted(genomes, key=itemgetter('Organism/Name'))

    # Semi-touch genomes file in case no genomes were selected, for instance when uploading external genomes
    open(genomes_file, mode='a').close()

    # Write IDs to file, with organism name as second column to make the project ID files more self explanatory.
    for genome in genomes:
        # Download files here, but ignore returned files: These can be retrieved from cache during extraction/translation
        download_genome_files(genome, genomes_file, require_ptt=require_ptt)

    # Post check after translation to see if more than one genome actually had some genomic contents
    with open(genomes_file) as read_handle:
        genome_ids = [line.split()[0] for line in read_handle]
        # If some genomes were skipped, ensure at least two genomes remain
        if len([gid for gid in genome_ids if gid.startswith('#')]):
            assert 2 <= len([gid for gid in genome_ids if not gid.startswith('#')]), \
                "Some genomes were skipped, leaving us with less than two genomes to operate on; " \
                "Inspect messages in Project ID list and reevaluate genome selection"

    # Exit after a comforting log message
    logging.info("Produced: \n%s", genomes_file)

if __name__ == '__main__':
    main(sys.argv[1:])
