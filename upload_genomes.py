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
"""Module for the select taxa step."""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from divergence import parse_options, create_archive_of_files
import logging as log
import os
import sys
import tempfile

def format_fasta_genome_headers(label, nucl_fasta_file):
    """Format an individual nucleotide fasta file containing coding regions so the headers match expected patterns."""
    label = label.replace(' ', '_')

    #Determine output file name
    filename = os.path.split(nucl_fasta_file)[1]
    formatted_fasta_file = tempfile.mkstemp(suffix = '.ffn', prefix = '{0}.{1}.formatted_'.format(label, filename))[1]
    with open(formatted_fasta_file, mode = 'w') as write_handle:
        for index, nucl_seqrecord in enumerate(SeqIO.parse(nucl_fasta_file, 'fasta'), 1):
            #Try to retain user specified protein identifiers, hoping they stuck to this guideline:
            #http://www.ncbi.nlm.nih.gov/books/NBK7183/?rendertype=table&id=ch_demo.T5
            #But do throw in a little counter to make sure the generated IDs are actually unique
            protein_id = '{0:06}_{1}'.format(index, nucl_seqrecord.id.replace('|', '_'))

            #Alternatively create a new unique protein id using the enumerated number, first twelve bases and length
            #protein_id = 'protein_{0}:{1}...({2})'.format(index, str(nucl_sequence)[:12], len(nucl_sequence))

            #Remove gap codons from input nucleotide fasta file in complete codons only
            nucl_sequence_str = str(nucl_seqrecord.seq).replace('---', '')
            nucl_sequence = Seq(nucl_sequence_str)

            #Write out fasta. Header format as requested: >project_id|genbank_ac|protein_id|cog|source 
            header = '{0}|{1}|{2}|{3}|{4}'.format(label, filename, protein_id, None, 'upload')

            #Create protein sequence record and write it to file
            formatted_seqrecord = SeqRecord(nucl_sequence, id = header, description = '')
            SeqIO.write(formatted_seqrecord, write_handle, 'fasta')
    return formatted_fasta_file

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: select_taxa.py
--external-genomes=    comma-separated list of label:nucleotide fasta file pairs of externally supplied genomes.
    label:FILE,...     labels should be unique as genomes will be identified by this label in further output files
--external-zip=FILE    destination path for archive of user provided external genomes containing formatted nucleotide fasta files
"""
    options = ['external-genomes', 'external-zip']
    external_genomes, external_zip = parse_options(usage, options, args)

    #External genomes are nucleotide fasta files uploaded by the user of which we will reformat the header
    external_fasta_files = {}

    #Handle externally uploaded genomes
    #Sample line: label1:file1,label2:file2, #Note trailing the trailing , that's a Galaxy artifact we'll ignore
    for label, filename in (label_file.split(':') for label_file in external_genomes.split(',') if label_file):
        log.info('Formatting external genome labeled %s at %s', label, filename)
        formatted_file = format_fasta_genome_headers(label, filename)
        external_fasta_files[label] = formatted_file

    #Copy formatted external genome files to archive that will be output as well
    create_archive_of_files(external_zip, external_fasta_files.values())

    #Remove temporary formatted files
    for formatted_file in external_fasta_files.values():
        os.remove(formatted_file)

    #Exit after a comforting log message
    log.info("Produced: \n%s", external_zip)

if __name__ == '__main__':
    main(sys.argv[1:])
