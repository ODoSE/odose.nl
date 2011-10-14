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
"""Module to create a crosstable between orthologs & genomes showing gene IDs at intersections."""

from Bio import SeqIO
from divergence import parse_options, extract_archive_of_files
from operator import itemgetter
import logging
import os.path
import shutil
import sys
import tempfile

def create_crosstable(sico_files, target_crosstable):
    """Create crosstable with vertically the orthologs, horizontally the genomes, and gene IDs at intersections."""
    with open(target_crosstable, mode = 'w') as write_handle:
        #Create dictionaries mapping genomes to gene IDs per sico file 
        row_data = [(sico_file, dict(itemgetter(0, 2)(fasta_record.id.split('|'))
                         for fasta_record in SeqIO.parse(sico_file, 'fasta')))
                    for sico_file in sico_files]

        #Retreve unique genomes across all sico files, just to be safe 
        genomes = sorted(set(key for row in row_data for key in row[1].keys()))

        #Write out values to file
        write_handle.write('\t' + '\t'.join(genomes) + '\n')
        for sico_file, row in row_data:
            ortholog = os.path.split(sico_file)[1].split('.')[0]
            write_handle.write(ortholog + '\t')
            write_handle.write('\t'.join(row.get(genome, '') for genome in genomes))
            write_handle.write('\n')

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: crosstable_gene_ids.py
--sico-zip=FILE      archive of single copy orthologous (SICO) genes in separate files per ortholog
--crosstable=FILE    destination file path for crosstable between orthologs & genomes with gene IDs at intersections
"""
    options = ['sico-zip', 'crosstable']
    sizo_zip, target_crosstable = parse_options(usage, options, args)

    #Create tempdir
    run_dir = tempfile.mkdtemp(prefix = 'crosstable_')
    sico_files = extract_archive_of_files(sizo_zip, run_dir)

    #Create crosstable
    create_crosstable(sico_files, target_crosstable)

    #Remove extracted files to free disk space 
    shutil.rmtree(run_dir)

    #Exit after a comforting log message
    logging.info("Produced: \n%s", target_crosstable)

if __name__ == '__main__':
    main(sys.argv[1:])
