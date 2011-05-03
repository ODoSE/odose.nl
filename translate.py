#!/usr/bin/env python
"""Module to translate between DNA and protein level."""

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from divergence import create_directory, concatenate, create_archive_of_files
from divergence.select_taxa import download_genome_files, select_genomes
from multiprocessing import Pool
import getopt
import logging as log
import os
import sys

def _build_protein_to_cog_mapping(ptt_file):
    """Build a dictionary mapping PID to COG based on protein table file."""
    mapping = {}
    with open(ptt_file) as read_handle:
        #Line 1: Escherichia coli 536, complete genome - 1..4938920
        read_handle.readline().strip() #Intentionally ignored
        #Line 2: 4619 proteins
        line = read_handle.readline().strip()
        assert line.endswith(" proteins")
        #Line 3: Location    Strand    Length    PID    Gene    Synonym    Code    COG    Product
        line = read_handle.readline().strip()
        columns = line.split("\t")
        assert columns == ['Location', 'Strand', 'Length', 'PID', 'Gene', 'Synonym', 'Code', 'COG', 'Product']

        for line in read_handle:
            values = line.split('\t')
            assert len(values) == 9
            #Make format match format in genbank db_xref records
            pid = 'GI:' + values[3]
            assert pid not in mapping, 'Protein identifier was already assigned value: ' + mapping[pid]
            #Assign None if value of COG does not start with COG (such as when it's '-')
            cog = values[7]
            mapping[pid] = cog if cog.startswith('COG') else None
    return mapping

def translate_genomes(genomes):
    """Download genome files, extract genes and translate those to proteins, returning DNA and protein fasta files."""
    assert len(genomes), 'Some genomes should be selected'

    pool = Pool()
    futures = [(genome, pool.apply_async(download_genome_files, (genome,))) for genome in genomes]
    dna_aa_pairs = [_translate_genome(genome, gbk_ptt_pairs.get()) for genome, gbk_ptt_pairs in futures]
    #Unzip the above produced pairs using some magic, as per: http://docs.python.org/library/functions.html#zip
    dna_files, aa_files = zip(*dna_aa_pairs)
    return dna_files, aa_files

def _translate_genome(genome, tuples_of_gbk_and_ptt_files):
    """Translate all files for genome and concatenate them into single DNA and Protein fasta files."""
    assert tuples_of_gbk_and_ptt_files is not None, 'No genbank files were provided for genome: ' + genome

    refseq_id = genome['RefSeq project ID']
    out_dir = create_directory('translations/' + refseq_id)
    dna_files = []
    protein_files = []
    for gbk_file, ptt_file in tuples_of_gbk_and_ptt_files:
        dna_file, protein_file = _extract_gene_and_protein(out_dir, refseq_id, gbk_file, ptt_file)
        dna_files.append(dna_file)
        protein_files.append(protein_file)

    #Concatenate files into one
    dna_concatemer = os.path.join(out_dir, '{pid}.ffn'.format(pid = refseq_id))
    protein_concatemer = os.path.join(out_dir, '{pid}.faa'.format(pid = refseq_id))
    concatenate(dna_concatemer, dna_files)
    concatenate(protein_concatemer, protein_files)
    return dna_concatemer, protein_concatemer

def _extract_gene_and_protein(out_dir, refseq_id, genbank_file, ptt_file = None):
    """Translate genbank DNA to protein and return resulting fasta files."""
    #Determine filename for output file
    file_root = os.path.splitext(genbank_file)[0]
    aa_file = os.path.join(out_dir, os.path.split(file_root)[1] + '.faa')
    dna_file = os.path.join(out_dir, os.path.split(file_root)[1] + '.ffn')

    #Check if file already exists and has content, and if so return it
    if os.path.isfile(aa_file) and 0 < os.path.getsize(aa_file) and \
        os.path.isfile(dna_file) and 0 < os.path.getsize(dna_file):
        return dna_file, aa_file

    #Determine if we'll add COG annotation from the ptt file
    if ptt_file is None:
        log.warn('No .ptt file given for %s: COG assignment disabled', genbank_file)
        #Disable COG mapping by setting cog_dict to an empty dictionary
        cog_dict = {}
    else:
        #Retrieve PID to COG mapping from ptt file
        cog_dict = _build_protein_to_cog_mapping(ptt_file)

    #Open genbank_file & convert it using BioPython
    log.info('Translating %s', genbank_file)
    with open(genbank_file) as file_handle:
        with open(aa_file, mode = 'w') as aa_writer:
            with open(dna_file, mode = 'w') as dna_writer:
                for gb_record in SeqIO.parse(file_handle, 'genbank'):#Bio.GenBank.Record
                    plasmid = False
                    for gb_feature in gb_record.features:#Bio.SeqFeature
                        if gb_feature.type == 'source' and 'plasmid' in gb_feature.qualifiers:
                            plasmid = True
                        #Skip any non coding sequence features or pseudo (non-functional version) CDS
                        if gb_feature.type == 'CDS' and not 'pseudo' in gb_feature.qualifiers:
                            _extract_and_translate_cds(cog_dict, aa_writer, dna_writer, refseq_id, gb_record, gb_feature, plasmid)

    assert os.path.isfile(aa_file) and 0 < os.path.getsize(aa_file), 'File should exist and have some content now'
    assert os.path.isfile(dna_file) and 0 < os.path.getsize(dna_file), 'File should exist and have some content now'
    return dna_file, aa_file

def _extract_and_translate_cds(cog_mapping, aa_writer, dna_writer, refseq_id, gb_record, gb_feature, plasmid):
    """Extract DNA for Coding sequences, translate using GBK translation table and return dna & protein fasta files."""

    #Protein identifier is a property of the genbank feature
    protein_id = gb_feature.qualifiers['protein_id'][0]

    #Original sequence retrieved through BioPython 1.53+'s internal method
    extracted_seq = gb_feature.extract(gb_record.seq) #Bio.Seq.Seq

    #Translation table is a property of the genbank feature
    transl_table = gb_feature.qualifiers['transl_table'][0]

    #Translate entire sequence as coding sequence using above translation table
    #Additional CodonTables are optionally available from Bio.Data.CodonTable
    try:
        protein_seq = extracted_seq.translate(table = transl_table, cds = True)
    except TranslationError as err:
        if 'transl_except' in gb_feature.qualifiers:
            #Fall back on GenBank translation whenever a transl_except record is found
            protein_seq = gb_feature.qualifiers['translation'][0]
        else:
            log.error(gb_feature)
            log.error('%s: Error in translating %s\n%s', gb_record.id, protein_id, extracted_seq)
            raise err #Log some debug information before reraising error

    #Determine COG by looking it up based on protein identifier
    cog = None
    for xref in gb_feature.qualifiers['db_xref']:
        if xref in cog_mapping:
            cog = cog_mapping[xref]
            break

    #Strip stopcodon from DNA sequence here, so it's removed by the time we align & trim the sequences
    last_codon = str(extracted_seq[-3:])
    assert last_codon in ['TAA', 'TAG', 'TGA'], 'Expected stopcodon, after using cds=True above, but was ' + last_codon
    extracted_seq = extracted_seq[:-3]

    #Write out fasta. Header format as requested: >refseq_id|genbank_ac|protein_id|cog|source (either core or plasmid)
    source = 'plasmid' if plasmid else 'core'
    header = '{0}|{1}|{2}|{3}|{4}'.format(refseq_id, gb_record.id, protein_id, cog, source)
    _write_fasta(aa_writer, header, str(protein_seq))
    _write_fasta(dna_writer, header, str(extracted_seq))

def _write_fasta(write_handle, headerline, sequence):
    """Wrapper to write fasta header and sequence to file in correct format with new lines and limited line length."""
    #Write fasta to output file
    write_handle.write('>' + headerline + '\n')
    #Print rows with a length of 70 characters  
    while 70 < len(sequence):
        write_handle.write(sequence[:70] + '\n')
        sequence = sequence[70:]
    #Write any remainder
    if sequence:
        write_handle.write(sequence + '\n')

def main(args):
    """Main function called when run from command line or as part of pipeline."""

    def _parse_options(args):
        """Use getopt to parse command line argument options"""

        def _usage():
            """Print _usage information"""
            print """
Usage: translate.py 
--genomes=FILE        file with refseq id from complete genomes table on each line 
--dna-zip=FILE        destination file path for zip archive of extracted DNA files
--protein-zip=FILE    destination file path for zip archive of translated protein files
            """

        options = ['genomes', 'dna-zip', 'protein-zip']
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

    genomes_file, dna_zipfile, protein_zipfile = _parse_options(args)

    #Read genomes ids from genomes_file, each on their own line
    with open(genomes_file, mode = 'r') as read_handle:
        genome_ids = [line.strip() for line in read_handle]

    #Retrieve full genome dictionaries from complete genomes table in select_taxa
    genomes = select_genomes(genome_ids)

    #Actually translate the genomes to produced a set of files for both  dna files & protein files
    dna_files, protein_files = translate_genomes(genomes)

    #Write the produced files to command line argument filenames
    create_archive_of_files(dna_zipfile, dna_files)
    create_archive_of_files(protein_zipfile, protein_files)

    #Exit after a comforting log message
    log.info("Produced: \n%s &\n%s", dna_zipfile, protein_zipfile)

if __name__ == '__main__':
    main(sys.argv[1:])
