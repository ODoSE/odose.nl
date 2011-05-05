#!/usr/bin/env python
"""Module to run Phylogenetic Analysis by Maximum Likelihood (PAML)."""

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from divergence import create_directory, extract_archive_of_files, create_archive_of_files
from divergence.select_taxa import select_genomes_from_file
from subprocess import check_call, STDOUT
import getopt
import logging as log
import os.path
import sys
import tempfile

def run_paml(genomes_a, genomes_b, sico_files):
    """Run PAML for representatives of clades A and B in each of the SICO files, to calculate dN/dS."""
    #PAML runs inside a temporary folder, to prevent interference with simultaneous runs 
    paml_dir = tempfile.mkdtemp(prefix = 'paml_run_')

    #Pick the first genomes as representatives for each clade
    representative_a = genomes_a[0]['RefSeq project ID']
    representative_b = genomes_b[0]['RefSeq project ID']

    log.info('Running PAML for {0} aligned and trimmed SICOs'.format(len(sico_files)))
    return [_run_yn00(paml_dir, representative_a, representative_b, sico_file) for sico_file in sico_files]

YN00 = '/projects/divergence/software/paml44/bin/yn00'

def _run_yn00(paml_dir, repr_id_a, repr_id_b, sico_file):
    """Run yn00 from PAML for selected sequence records from sico_file, returning main nexus output file."""
    #Find sequences from above chosen clade representatives in each SICO file
    seqr_a = None
    seqr_b = None
    for seqrecord in SeqIO.parse(sico_file, 'fasta', alphabet = generic_dna):
        refseq_id = seqrecord.id.split('|')[0]
        if refseq_id == repr_id_a:
            seqr_a = seqrecord
        elif refseq_id == repr_id_b:
            seqr_b = seqrecord
        if seqr_a and seqr_b:
            break
    else:
        assert False, '{0}: Both seqr_a & seqr_b should have been set\n{1}\n{2}'.format(sico_file, seqr_a, seqr_b)

    #Write the representative sequence records out to file in yn00 compatible format
    sico = os.path.splitext(os.path.split(sico_file)[1])[0]
    yn00_dir = create_directory(sico, inside_dir = paml_dir)
    nexus_file = os.path.join(yn00_dir, sico + '.nexus')
    _write_nexus_file(seqr_a, seqr_b, nexus_file)

    #Generate yn00 configuration file
    output_file = os.path.join(yn00_dir, 'yn' + sico)
    config_file = os.path.join(yn00_dir, 'yn00.ctl')
    _write_config_file(nexus_file, output_file, config_file)

    #Run yn00
    command = [YN00, os.path.split(config_file)[1]]
    check_call(command, cwd = yn00_dir, stdout = open('/dev/null', mode = 'w'), stderr = STDOUT)

    assert os.path.isfile(output_file) and os.path.getsize(output_file), 'Expected some content in ' + output_file
    return output_file

def _write_nexus_file(seqr_a, seqr_b, nexus_file):
    """Write representative sequences out to a file in the PAML compatible nexus format."""
    nexus_contents = """
#NEXUS
begin data; 
   dimensions ntax={ntax} nchar={nchar}; 
   format datatype=dna missing=? gap=-; 
matrix 
clade_a  {seqa}
clade_b  {seqb}
;
end;
""".format(ntax = 2, nchar = len(seqr_a), seqa = str(seqr_a.seq), seqb = str(seqr_b.seq))
    with open(nexus_file, mode = 'w') as write_handle:
        #SeqIO.write([seqr_a, seqr_b], write_handle, 'nexus')
        write_handle.write(nexus_contents)

def _write_config_file(nexus_file, output_file, config_file):
    """Write a yn00 configuration file using relative paths to the nexus file and output file."""
    config_contents = """
   seqfile = {0} * sequence data file name
   outfile = {1} * main result file
   verbose = 0  * 1: detailed output (list sequences), 0: concise output

     icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below

 weighting = 0  * weighting pathways between codons (0/1)?
commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)? 
*    ndata = 1

* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
""".format(os.path.split(nexus_file)[1], os.path.split(output_file)[1])
    with open(config_file, mode = 'w') as write_handle:
        write_handle.write(config_contents)

def main(args):
    """Main function called when run from command line or as part of pipeline."""

    def _parse_options(args):
        """Use getopt to parse command line argument options"""

        def _usage():
            """Print _usage information"""
            print """
Usage: run_paml.py 
--genomes-a=FILE    file with RefSeq id from complete genomes table on each line for clade A
--genomes-b=FILE    file with RefSeq id from complete genomes table on each line for clade B
--sico-zip=FILE     archive of aligned & trimmed single copy orthologous (SICO) genes
--paml-zip=FILE     destination file path for archive of PAML output per SICO gene
"""

        options = ['genomes-a', 'genomes-b', 'sico-zip', 'paml-zip']
        try:
            #postfix '=' to indicate options require an argument
            long_options = [opt + '=' for opt in options]
            tuples = getopt.getopt(args, '', long_options)[0]
            arguments = dict((opt[2:], value) for opt, value in tuples)
        except getopt.GetoptError as err:
            print str(err)
            _usage()
            sys.exit(1)

        #Ensure all arguments were provided
        for opt in options:
            if opt not in arguments:
                print 'Mandatory argument {0} not provided'.format(opt)
                _usage()
                sys.exit(1)

        #Retrieve & return file paths from dictionary
        return [arguments[option] for option in options]

    genome_a_ids_file, genome_b_ids_file, sico_zip, paml_zip = _parse_options(args)

    #Parse file containing RefSeq project IDs & retrieve associated genome dictionaries from complete genomes table
    genomes_a = select_genomes_from_file(genome_a_ids_file)
    genomes_b = select_genomes_from_file(genome_b_ids_file)

    #Extract files from zip archive
    temp_dir = tempfile.mkdtemp()
    sico_files = extract_archive_of_files(sico_zip, temp_dir)

    #Actually run cleanup
    paml_files = run_paml(genomes_a, genomes_b, sico_files)

    #Write the produced files to command line argument filenames
    create_archive_of_files(paml_zip, paml_files)

    #Exit after a comforting log message
    log.info("Produced: \n%s", paml_zip)

if __name__ == '__main__':
    main(sys.argv[1:])
