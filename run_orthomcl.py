#!/usr/bin/env python
"""Module to run orthoMCL. Steps in this module reflect the steps in the UserGuide.txt bundled with OrthoMCL."""

from Bio import SeqIO
from divergence import create_directory, extract_archive_of_files, parse_options
from divergence.orthomcl_database import create_database, get_configuration_file, delete_database
from divergence.reciprocal_blast_lsgp import reciprocal_blast  # Remove _lsgp suffix to using local BLAST instead
from divergence.translate import translate_fasta_coding_regions
from divergence.upload_genomes import format_fasta_genome_headers
from divergence.versions import MCL, ORTHOMCL_INSTALL_SCHEMA, ORTHOMCL_ADJUST_FASTA, ORTHOMCL_FILTER_FASTA, \
    ORTHOMCL_BLAST_PARSER, ORTHOMCL_LOAD_BLAST, ORTHOMCL_PAIRS, ORTHOMCL_DUMP_PAIRS_FILES
from subprocess import Popen, PIPE, CalledProcessError, check_call, STDOUT
import logging as log
import multiprocessing
import os
import shutil
import sys
import tempfile


__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


def run_orthomcl(proteome_files, poor_protein_length, evalue_exponent, target_poor_proteins_file, target_groups_file):
    """Run all the steps in the orthomcl pipeline, starting with a set of proteomes and ending up with groups.txt."""
    #Delete orthomcl directory to prevent lingering files from previous runs to influence new runs
    run_dir = tempfile.mkdtemp(prefix='orthomcl_run_')

    #Steps leading up to and performing the reciprocal blast, as well as minor post processing
    adjusted_fasta_dir, fasta_files = _step5_orthomcl_adjust_fasta(run_dir, proteome_files)
    good, poor = _step6_orthomcl_filter_fasta(run_dir, adjusted_fasta_dir, min_length=poor_protein_length)
    allvsall = _step7_blast_all_vs_all(good, fasta_files)
    similar_sequences = _step8_orthomcl_blast_parser(run_dir, allvsall, adjusted_fasta_dir)
    #Clean up all vs all blast results file
    os.remove(allvsall)

    #Create new database and install database schema in it, so individual runs do not interfere with each other
    dbname = create_database()
    config_file = get_configuration_file(run_dir, dbname, evalue_exponent)
    _step4_orthomcl_install_schema(run_dir, config_file)

    #Steps that occur in database, and thus do little to produce output files
    _step9_orthomcl_load_blast(similar_sequences, config_file)
    _step10_orthomcl_pairs(run_dir, config_file)
    mcl_input = _step11_orthomcl_dump_pairs(run_dir, config_file)[0]

    #Trash database now that we're done with it
    delete_database(dbname)

    #MCL related steps: run MCL on mcl_input resulting in the groups.txt file
    groups = _step12_mcl(run_dir, mcl_input)

    #Move poor proteins file & groups file outside run_dir ahead of removing run_dir
    shutil.move(poor, target_poor_proteins_file)
    shutil.move(groups, target_groups_file)
    #Remove run_dir to free disk space
    shutil.rmtree(run_dir)

    return target_groups_file, target_poor_proteins_file


def _step4_orthomcl_install_schema(run_dir, config_file):
    """Create OrthoMCL schema in an Oracle or Mysql database.

    usage: orthomclInstallSchema config_file sql_log_file

    where:
      config_file : orthomcl configuration file
      sql_log_file : optional log of sql executed

    EXAMPLE: orthomclSoftware/bin/orthomclInstallSchema my_orthomcl_dir/orthomcl.config my_orthomcl_dir/install_schema.log

    NOTE: the database login in the config file must have update/insert/truncate privileges on the tables specified in the config file.
    """
    sql_log_file = tempfile.mkstemp('.log', 'orthomclInstallSchema_run_', dir=run_dir)[1]
    command = [ORTHOMCL_INSTALL_SCHEMA, config_file, sql_log_file]
    log.info('Executing: %s', ' '.join(command))
    check_call(command)
    return sql_log_file


def _step5_orthomcl_adjust_fasta(run_dir, proteome_files, id_field=3):
    """Create an OrthoMCL compliant .fasta file, by adjusting definition lines.

    Usage:
      orthomclAdjustFasta taxon_code fasta_file id_field

    where:
      taxon_code:  a three or four letter unique abbreviation for the taxon
      fasta_file:  the input fasta file per proteome
      id_field:    a number indicating what field in the definition line contains
                   the protein ID.  Fields are separated by either ' ' or '|'. Any
                   spaces immediately following the '>' are ignored.  The first
                   field is 1. For example, in the following definition line, the
                   ID (AP_000668.1) is in field 4:  >gi|89106888|ref|AP_000668.1|

    Input file requirements:
      (1) .fasta format
      (2) a unique id is provided for each sequence, and is in the field specified
          by id_field

    Output file format:
      (1) .fasta format
      (2) definition line is of the form:
             >taxoncode|unique_protein_id

    The output file is named taxoncode.fasta

    Note: if your input files do not meet the requirements, you can do some simple perl or awk processing of them to
    create the required input files to this program, or the required output files.  This program is provided as a
    convenience, but OrthoMCL users are expected to have the scripting skills to provide compliant .fasta files.

    EXAMPLE: orthomclSoftware/bin/orthomclAdjustFasta hsa Homo_sapiens.NCBI36.53.pep.all.fa 1
    """
    #Create directory to hold compliant fasta
    adjusted_fasta_dir = create_directory('compliant_fasta', inside_dir=run_dir)
    adjusted_fasta_files = []
    for proteome_file in proteome_files:
        #Use first part of header of first entry as taxon code
        record_iter = SeqIO.parse(proteome_file, 'fasta')
        taxon_code = record_iter.next().id.split('|')[0]  # pylint: disable=E1101

        #Call orhtomclAdjustFasta
        command = [ORTHOMCL_ADJUST_FASTA, taxon_code, proteome_file, str(id_field)]
        log.info('Executing: %s', ' '.join(command))
        check_call(command)
        #Move resulting fasta file to compliantFasta directory
        adjusted_fasta_file = taxon_code + '.fasta'
        fasta_file_destination = os.path.join(adjusted_fasta_dir, adjusted_fasta_file)
        shutil.move(adjusted_fasta_file, fasta_file_destination)
        adjusted_fasta_files.append(fasta_file_destination)
    #Return path to directory containing compliantFasta
    return adjusted_fasta_dir, adjusted_fasta_files


def _step6_orthomcl_filter_fasta(run_dir, input_dir, min_length=10, max_percent_stop=20):
    """Create goodProteins.fasta containing all good proteins and rejectProteins.fasta containing all rejects. Input is
    a directory containing a set of compliant input .fasta files (as produced by orthomclAdjustFasta).

    Usage:
      orthomclFilterFasta input_dir min_length max_percent_stops

    where:
      input_dir:           a directory containing a set of .fasta files
      min_length:          minimum allowed length of proteins.  (suggested: 10)
      max_percent_stop:    maximum percent stop codons.  (suggested 20)

    The input requirements are:
      1) a compliantFasta/ directory which contains all and only the proteome .fasta files, one file per proteome.
      2) each .fasta file must have a name in the form 'xxxx.fasta' where xxxx is a three or four letter unique taxon
         code.  For example: hsa.fasta or eco.fasta
      3) each protein in those files must have a definition line in the following format:
         >xxxx|yyyyyy
         where xxxx is the three or four letter taxon code and yyyyyy is a sequence identifier unique within that taxon.

    Output:
      - my_orthomcl_dir/goodProteins.fasta
      - my_orthomcl_dir/poorProteins.fasta
      - report of suspicious proteomes (> 10% poor proteins)

    EXAMPLE: orthomclSoftware/bin/orthomclFilterFasta my_orthomcl_dir/compliantFasta 10 20
    """
    #Run orthomclFilterFasta
    out_dir = create_directory('filtered_fasta', inside_dir=run_dir)
    report = os.path.join(out_dir, 'filter_report.log')
    with open(report, mode='w') as report_file:
        command = [ORTHOMCL_FILTER_FASTA, input_dir, str(min_length), str(max_percent_stop)]
        log.info('Executing: %s', ' '.join(command))
        check_call(command, stdout=report_file, stderr=STDOUT)

    #Move output files to out directory
    good = os.path.join(out_dir, 'good_proteins.fasta')
    poor = os.path.join(out_dir, 'poor_proteins.fasta')
    shutil.move('goodProteins.fasta', good)
    shutil.move('poorProteins.fasta', poor)

    #Ensure neither of the proteomes is suspicious according to min_length & max_percent_stop
    with open(report) as report_file:
        if 'Proteomes with > 10% poor proteins:' in report_file.read():  # OrthoMCL does NOT add actual min_length value
            msg = 'OrthomclFilterFasta found suspicious proteomes based on values for length'
            log.error(msg)
            assert False, msg

    #Warn the user about the poor proteins found here, if they were found at all
    poor_records = list(SeqIO.parse(poor, 'fasta'))
    if poor_records:
        log.warn('%i poor sequence records identified by orthomclFilterFasta:', len(poor_records))
        for seqr in poor_records:
            log.warn('>%s: %s', seqr.id, seqr.seq)

    #Assert good exists and has some content
    assert os.path.isfile(good) and 0 < os.path.getsize(good), good + ' should exist and have some content'

    #Only good and poor proteins
    return good, poor


def _step7_blast_all_vs_all(good_proteins_file, fasta_files):
    """Input:
      - goodProteins.fasta
    Output:
      - your_blast_results_in_tab_format

    You must run your own BLAST.  For large datasets you should consider gaining access to a compute cluster.

    We expect you to:
      - use NCBI BLAST
      - run with the -m 8 option to provide tab delimited output required by Step 8
      - for IMPORTANT DETAILS about other BLAST arguments, see:
        the OrthoMCL Algorithm Document (http://docs.google.com/Doc?id=dd996jxg_1gsqsp6)

    If you are a power user you can deviate from this, so long as you can ultimately provide output in exactly the format provided by NCBI BLAST using the -m 8 option, and expected by Step 8.

    If you are a super - power user you can deviate from that, and also skip Step 8.   But you must be able to provide the exact format file created by that step as expected by Step 9.  The tricky part is computing percent match.

    Time estimate: highly dependent on your data and hardware
    """
    #Handled by reciprocal blast module
    return reciprocal_blast(good_proteins_file, fasta_files)


def _step8_orthomcl_blast_parser(run_dir, blast_file, fasta_files_dir):
    """orthomclBlastParser blast_file fasta_files_dir

    where:
      blast_file:       BLAST output in m8 format.
      fasta_files_dir:  a directory of compliant fasta files as produced by
                        orthomclAdjustFasta

    m8 format has these columns:
      query_name, hitname, pcid, len, mismatches, ngaps, start('query'),
      end('query'), start('hit'), end('hit'), evalue, bits

    output:
      tab delimited text file, with one row per query-subject match. the columns are:
         query_id, subject_id, query_taxon, subject_taxon,
         evalue_mant, evalue_exp, percent_ident, percent_match

    (percent_match is computed by counting the number of bases or amino acids in the shorter sequence that are matched in any hsp, and dividing by the length of that shorter sequence)

    EXAMPLE: orthomclSoftware/bin/orthomclBlastParser my_blast_results my_orthomcl_dir/compliantFasta >> my_orthomcl_dir/similar_sequences.txt
    """
    #Run orthomclBlastParser
    command = [ORTHOMCL_BLAST_PARSER, blast_file, fasta_files_dir]
    log.info('Executing: %s', ' '.join(command))
    similar_sequences = os.path.join(run_dir, 'similar_sequences.tsv')
    with open(similar_sequences, mode='w') as stdout_file:
        #check_call(command, stdout = stdout_file, stderr = open('/dev/null', mode = 'w'))
        process = Popen(command, stdout=stdout_file, stderr=PIPE)
        retcode = process.wait()
        if retcode:
            stderr = process.communicate()[1]
            log.error(stderr)
            raise CalledProcessError(retcode, command)

    msg = 'Similar seqeunces files should now have some content'
    assert os.path.isfile(similar_sequences) and 0 < os.path.getsize(similar_sequences), msg

    return similar_sequences


#Steps 9, 10 and 11 all use the same relational database, which could cause problems with simultaneous runs
def _step9_orthomcl_load_blast(similar_seqs_file, config_file):
    """Load Blast results into an Oracle or Mysql database.

    usage: orthomclLoadBlast config_file similar_seqs_file

    where:
      config_file :       see below
      similar_seqs_file : output from orthomclParseBlast

    EXAMPLE: orthomclSoftware/bin/orthomclLoadBlast my_orthomcl_dir/orthomcl.config my_orthomcl_dir/similar_sequences.txt
    """
    #Run orthomclLoadBlast
    command = [ORTHOMCL_LOAD_BLAST, config_file, similar_seqs_file]
    log.info('Executing: %s', ' '.join(command))
    check_call(command)
    return


def _step10_orthomcl_pairs(run_dir, config_file):
    """Find pairs for OrthoMCL.

    usage: orthomclPairs config_file log_file cleanup=[yes|no|only|all] <startAfter=TAG>

    where:
      config_file : see below
      cleanup     : clean up temp tables?
                       yes = clean as we go;
                       no = don't clean as we go;
                       only=just clean, do nothing else;
                       all=just clean, plus clean InParalog, Ortholog and CoOrtholog tables.
      startAfter  : optionally start after a previously completed step. see below for TAGs

    Database Input:
      - SimilarSequences table containing all-v-all BLAST hits
      - InParalog, Ortholog, CoOrtholog tables - created but empty

    Database Output:
      - Populated InParalog, Ortholog and CoOrtholog tables

    Standard Error:
      - logging info

    NOTE: the database login in the config file must have update/insert/truncate privileges on the tables specified in the config file.

    EXAMPLE: orthomclSoftware/bin/orthomclPairs my_orthomcl_dir/orthomcl.config my_orthomcl_dir/orthomcl_pairs.log cleanup=no
    """
    #Run orthomclPairs
    pairs_log = os.path.join(run_dir, 'orthomclPairs.log')
    command = [ORTHOMCL_PAIRS, config_file, pairs_log, 'cleanup=no']
    log.info('Executing: %s', ' '.join(command))
    check_call(command)

    return pairs_log


def _step11_orthomcl_dump_pairs(run_dir, config_file):
    """Dump files from the database produced by the orthomclPairs program.

    usage: orthomclDumpPairsFiles config_file

    where:
      config_file : see below (you can use the same file given to orthomclPairs)

    Database Input:
      - InParalog, Ortholog, CoOrtholog tables - populated by orthomclPairs

    Output files:
      orthomclMclInput                       - file required by the mcl program
      pairs/                                 - dir holding relationship files
        potentialOrthologs.txt               - ortholog relationships
        potentialInparalogs.txt              - inparalog relationships
        potentialCoorthologs.txt             - coortholog relationships

    The pairs/ files contain the pairs found by the orthomclPairs tables, and their
    average normalized scores.  This is the same information as in the
    orthomclMclInput file, but segregated by relationship type.  These are
    candidate relationships (edges) that will subsequently be grouped (clustered)
    by the mcl program to form the OrthoMCL ortholog groups.  These files contain
    more sensitive and less selective relationships then the final ortholog groups.

    Standard Error:
      - logging info

    EXAMPLE: orthomclSoftware/bin/orthomclDumpPairsFile out_dir/orthomcl.config
    """
    #Run orthomclDumpPairsFile
    out_dir = create_directory('orthologs', inside_dir=run_dir)
    command = [ORTHOMCL_DUMP_PAIRS_FILES, config_file]
    log.info('Executing: %s', ' '.join(command))
    check_call(command, cwd=out_dir)

    #Desired destination output file paths
    mcl_dir = create_directory('mcl', inside_dir=run_dir)
    mclinput = os.path.join(mcl_dir, 'mclInput.tsv')
    orthologs = os.path.join(out_dir, 'potentialOrthologs.tsv')
    inparalogs = os.path.join(out_dir, 'potentialInparalogs.tsv')
    coorthologs = os.path.join(out_dir, 'potentialCoorthologs.tsv')

    #Move output files to desired destinations
    shutil.move(os.path.join(out_dir, 'mclInput'), mclinput)
    shutil.move(os.path.join(out_dir, 'pairs/orthologs.txt'), orthologs)
    shutil.move(os.path.join(out_dir, 'pairs/inparalogs.txt'), inparalogs)
    shutil.move(os.path.join(out_dir, 'pairs/coorthologs.txt'), coorthologs)

    #Assert mcl input file exists and has some content
    assert os.path.isfile(mclinput) and 0 < os.path.getsize(mclinput), mclinput + ' should exist and have some content'

    return mclinput, orthologs, inparalogs, coorthologs


def _step12_mcl(run_dir, mcl_input_file):
    """Markov Cluster Algorithm: http://www.micans.org/mcl/

    Input:
      - mclInput file
    Output:
      - mclOutput file

    mcl my_orthomcl_dir/mclInput --abc -I 1.5 -o my_orthomcl_dir/mclOutput
    """
    #Run mcl
    mcl_dir = create_directory('mcl', inside_dir=run_dir)
    mcl_output_file = os.path.join(mcl_dir, 'mclOutput.tsv')
    mcl_log = os.path.join(mcl_dir, 'mcl.log')
    with open(mcl_log, mode='w') as open_file:
        threads = str(multiprocessing.cpu_count())
        command = [MCL, mcl_input_file, '--abc', '-I', '1.5', '-o', mcl_output_file, '-te', threads]
        log.info('Executing: %s', ' '.join(command))
        check_call(command, stdout=open_file, stderr=STDOUT)
    return mcl_output_file


def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: run_orthomcl.py
--protein-zip=FILE           zip archive of translated protein files
--ortholog-limiter=FILE      nucleotide fasta file containing coding regions in individual records. this file will be
                             translated to protein and fed into orthomcl along with files in protein - zip to influence
                             the clustering of orthologs, and (optionally) later the extraction of orthologs [OPTIONAL]
--poor-protein-length=INT    filter poor proteins when smaller than poor-protein-length
--evalue-exponent=INT        filter OrthoMCL BLAST similarities with Expect value exponents greater than this value
--poor-proteins=FILE         destination file path for filtered poor proteins
--groups=FILE                destination file path for file listing groups of orthologous proteins
"""
    options = ['protein-zip', 'ortholog-limiter=?', 'poor-protein-length', 'evalue-exponent', 'poor-proteins', 'groups']
    protein_zipfile, limiter_file, poor_protein_length, evalue_exponent, target_poor_proteins, target_groups_path = \
        parse_options(usage, options, args)

    #Extract files from zip archive
    temp_dir = tempfile.mkdtemp(prefix='orthomcl_proteins_')
    proteome_files = extract_archive_of_files(protein_zipfile, temp_dir)

    #If limiter file is defined, add it to the set op protein files
    if limiter_file:
        #First format nucleotide fasta file to contain the correct fasta headers
        formatted_fasta_file = format_fasta_genome_headers('limiter', limiter_file)
        #Then translate it from nucleotide to protein
        translated_limiter = translate_fasta_coding_regions(formatted_fasta_file)
        #Then append it to the list op proteome files
        proteome_files.append(translated_limiter)

    #Actually run orthomcl
    run_orthomcl(proteome_files, poor_protein_length, evalue_exponent, target_poor_proteins, target_groups_path)

    #Remove unused files to free disk space
    shutil.rmtree(temp_dir)

    #Exit after a comforting log message
    log.info("Produced: \n%s\n%s", target_poor_proteins, target_groups_path)

if __name__ == '__main__':
    main(sys.argv[1:])
