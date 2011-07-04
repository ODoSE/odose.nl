#!/usr/bin/env python
"""Module to translate between DNA and protein level."""

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Data.CodonTable import TranslationError
from divergence import create_directory, concatenate, create_archive_of_files, parse_options
from divergence.select_taxa import download_genome_files, select_genomes_by_ids
from multiprocessing import Pool
from operator import itemgetter
import logging as log
import os
import shutil
import sys
import tempfile

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

    #Use a pool to download files in the background while translating
    pool = Pool()
    futures = [pool.apply_async(download_genome_files, (genome,)) for genome in genomes]
    dna_aa_pairs = [_translate_genome(gbk_ptt_pairs.get()) for gbk_ptt_pairs in futures if gbk_ptt_pairs.get() != None]

    #Extract DNA & Protein files separately from dna_aa_pairs
    dna_files = [pair[0] for pair in dna_aa_pairs]
    aa_files = [pair[1] for pair in dna_aa_pairs]

    return dna_files, aa_files

def _translate_genome(tuples_of_gbk_and_ptt_files):
    """Translate all files for genome and concatenate them into single DNA and Protein fasta files."""
    assert tuples_of_gbk_and_ptt_files is not None, 'No genbank files were provided'

    project_id = tuples_of_gbk_and_ptt_files[0][0]
    out_dir = create_directory('translations/' + project_id)
    dna_files = []
    protein_files = []
    for project_id, gbk_file, ptt_file in tuples_of_gbk_and_ptt_files:
        dna_file, protein_file = _extract_gene_and_protein(out_dir, project_id, gbk_file, ptt_file)
        dna_files.append(dna_file)
        protein_files.append(protein_file)

    #Concatenate files into one
    dna_concatemer = os.path.join(out_dir, '{pid}.ffn'.format(pid = project_id))
    protein_concatemer = os.path.join(out_dir, '{pid}.faa'.format(pid = project_id))
    concatenate(dna_concatemer, dna_files)
    concatenate(protein_concatemer, protein_files)
    return dna_concatemer, protein_concatemer

def _extract_gene_and_protein(out_dir, project_id, genbank_file, ptt_file = None):
    """Translate genbank DNA to protein and return resulting fasta files."""
    #Determine filenames for temporary and cache destination output files
    file_root = os.path.splitext(genbank_file)[0]
    dna_file_dest = os.path.join(out_dir, os.path.split(file_root)[1] + '.ffn')
    aa_file_dest = os.path.join(out_dir, os.path.split(file_root)[1] + '.faa')

    #Check if destination output files already exist and have content, and if so return them
    if (os.path.isfile(aa_file_dest) and 0 < os.path.getsize(aa_file_dest) and
        os.path.isfile(dna_file_dest) and 0 < os.path.getsize(dna_file_dest)):
        return dna_file_dest, aa_file_dest

    #Use temporary files as write handles when translating, so we can not pollute cache with incomplete files
    dna_tmp = tempfile.mkstemp(suffix = '.ffn', prefix = 'translate_')[1]
    aa_tmp = tempfile.mkstemp(suffix = '.faa', prefix = 'translate_')[1]

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
    with open(aa_tmp, mode = 'w') as aa_wrtr:
        with open(dna_tmp, mode = 'w') as dna_wrtr:
            #Bio.GenBank.Record
            gb_recrd = SeqIO.read(genbank_file, 'genbank')

            #Determine whether the source of this genbank file is core genome or plasmid
            plasmid = False
            for gb_featr in gb_recrd.features:#Bio.SeqFeature
                if gb_featr.type == 'source' and 'plasmid' in gb_featr.qualifiers:
                    plasmid = True
                    break

            #Select only the coding sequences from all feature records
            coding_features = [gb_featr for gb_featr in gb_recrd.features
                               #Skip any non coding sequence features or pseudo (non-functional version) CDS
                               if gb_featr.type == 'CDS' and not 'pseudo' in gb_featr.qualifiers]

            #If there are no coding features, report this back to the user with a clear message rather than empty file
            if 0 == len(coding_features):
                log.error('No coding sequences found in genbank file %s', genbank_file)
                #TODO What are the consequences & how can the user recover from this?

            #Translate all coding sequences
            for cds_featr in coding_features:
                _extract_and_translate_cds(cog_dict, aa_wrtr, dna_wrtr, project_id, gb_recrd, cds_featr, plasmid)

    assert os.path.isfile(dna_tmp) and 0 < os.path.getsize(dna_tmp), dna_tmp + ' should exist and have some content'
    assert os.path.isfile(aa_tmp) and 0 < os.path.getsize(aa_tmp), aa_tmp + ' should exist and have some content'

    #Move completed files to cache location only just now, so incomplete files do not pollute cache when raising errors
    shutil.move(dna_tmp, dna_file_dest)
    shutil.move(aa_tmp, aa_file_dest)

    return dna_file_dest, aa_file_dest

def _extract_and_translate_cds(cog_mapping, aa_writer, dna_writer, project_id, gb_record, gb_feature, plasmid):
    """Extract DNA for Coding sequences, translate using GBK translation table and return dna & protein fasta files."""

    #Protein identifier is a property of the genbank feature
    protein_id = gb_feature.qualifiers['protein_id'][0]

    #Original sequence retrieved through BioPython 1.53+'s internal method
    extracted_seq = gb_feature.extract(gb_record.seq) #Bio.Seq.Seq

    if len(extracted_seq) % 3:
        #Skip CDS feature if length of extracted_seq is not a multiple of three
        log.warn('Length of extracted coding sequence %s not a multiple of 3: %i\n%s',
                  protein_id, len(extracted_seq), extracted_seq)
        return

    #Translation table is a property of the genbank feature
    transl_table = gb_feature.qualifiers['transl_table'][0]
    codon_table = CodonTable.unambiguous_dna_by_id.get(int(transl_table))

    #Set flag only when this CDS ends in a stop codon, so we can strip it off later, but do not strip non-stop-codons 
    cds_has_stopcodon = str(extracted_seq[-3:]) in codon_table.stop_codons

    #Translate entire sequence as coding sequence using above translation table
    #Additional CodonTables are optionally available from Bio.Data.CodonTable
    try:
        #If transl_except is present, an alternative translation is offered in the GenBank record
        if 'transl_except' in gb_feature.qualifiers:
            #Fall back on GenBank translation whenever a transl_except record is found
            protein_seq = gb_feature.qualifiers['translation'][0]
        else:
            protein_seq = extracted_seq.translate(table = transl_table, cds = True)
    except TranslationError as err:
        #Try to recover from the following errors
        if ('First codon ' in str(err) and ' is not a start codon' in str(err)) or \
           ('Final codon ' in str(err) and ' is not a stop codon' in str(err)):
            #Retrieve GenBank provided reference translation
            provided_seq = gb_feature.qualifiers['translation'][0]
            #Calculate lengths for DNA & protein sequences
            dna_length = len(extracted_seq)
            aa_length = len(provided_seq)

            #Occasionally an incomplete amino end coding sequence is found, such as in: YP_004377721.1 
            log.warning('Incomplete CDS found for %s from %s:\n%s', protein_id, gb_record.id, extracted_seq)

            #Calculate lengths for DNA & protein sequences
            if cds_has_stopcodon:
                #Correct check for the presence of a final stop codon
                dna_length -= 3

            #Assert lengths match, when stopcodon length is deducted
            if dna_length != 3 * aa_length:
                log.error('Triple protein length %i did not match DNA length %i', 3 * aa_length, dna_length)
                #Skip this CDS feature
                return

            protein_seq = provided_seq
            log.warning('Reverted to GenBank provided translation:\n%s', protein_seq)
        else:
            #Log some debug information before reraising error
            log.error('Error in translating %s from %s using:\n%s', protein_id, gb_record.id, extracted_seq)
            log.error(gb_feature)
            log.error(err)
            #Skip this CDS feature, but do continue translating future CDS within the same GenBank record 
            return

    #Determine COG by looking it up based on protein identifier
    cog = None
    for xref in gb_feature.qualifiers['db_xref']:
        if xref in cog_mapping:
            cog = cog_mapping[xref]
            break

    #Strip stopcodon from DNA sequence here, so it's removed by the time we align & trim the sequences
    if cds_has_stopcodon:
        extracted_seq = extracted_seq[:-3]

    #Write out fasta. Header format as requested: >project_id|genbank_ac|protein_id|cog|source (either core or plasmid)
    source = 'plasmid' if plasmid else 'core'
    header = '{0}|{1}|{2}|{3}|{4}'.format(project_id, gb_record.id, protein_id, cog, source)
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
    usage = """
Usage: translate.py 
--genomes=FILE        file with genbank id from complete genomes table on each line 
--dna-zip=FILE        destination file path for zip archive of extracted DNA files
--protein-zip=FILE    destination file path for zip archive of translated protein files
"""
    options = ['genomes', 'dna-zip', 'protein-zip']
    genome_ids_file, dna_zipfile, protein_zipfile = parse_options(usage, options, args)

    #Read GenBank Project IDs from genomes_file, each on their own line
    with open(genome_ids_file) as read_handle:
        genome_ids = [line.split()[0] for line in read_handle]

    #Retrieve associated genome dictionaries from complete genomes table
    genomes = select_genomes_by_ids(genome_ids).values()
    genomes = sorted(genomes, key = itemgetter('Organism Name'))

    #Actually translate the genomes to produced a set of files for both  dna files & protein files
    dna_files, protein_files = translate_genomes(genomes)

    #Write the produced files to command line argument filenames
    create_archive_of_files(dna_zipfile, dna_files)
    create_archive_of_files(protein_zipfile, protein_files)

    #Do not clean up extracted DNA files or Protein translations: Keep them as cache

    #Exit after a comforting log message
    log.info("Produced: \n%s &\n%s", dna_zipfile, protein_zipfile)

if __name__ == '__main__':
    main(sys.argv[1:])
