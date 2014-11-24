#!/usr/bin/env python
"""Module to translate between DNA and protein level."""

from Bio import SeqIO
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Data import CodonTable
from Bio.Data.CodonTable import TranslationError
from Bio.SeqRecord import SeqRecord
from divergence import create_directory, concatenate, create_archive_of_files, parse_options, \
    extract_archive_of_files, CODON_TABLE_ID
from divergence.select_taxa import select_genomes_by_ids
from divergence.download_taxa_mrs import download_genome_files, download_plasmid_files
from multiprocessing import Pool
from operator import itemgetter
import logging as log
import os
import re
import shutil
import sys
import tempfile

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"

# Using the standard NCBI Bacterial, Archaeal and Plant Plastid Code translation table (11).
BACTERIAL_CODON_TABLE = CodonTable.unambiguous_dna_by_id.get(CODON_TABLE_ID)


def _append_external_genomes(external_fasta_files, genomes_file):
    """Read out user provided labels and original filenames for uploaded genomes and append them to genome IDs file."""
    for fasta_file in external_fasta_files:
        # Read user supplied species label & originating filename
        with open(fasta_file) as read_handle:
            # Read first line, which should be header of first record, starting with a skipped >
            header = read_handle.readline()[1:]

            # Header format as requested: >project_id|genbank_ac|protein_id|cog|source
            label, origin = header.split('|')[0:2]

            with open(genomes_file, mode='a') as append_handle:
                # We'll use this 'external genome' source to skip externally derived files when downloading & translating
                append_handle.write('{0}\t{1}\texternal genome\n'.format(label, origin))


def _map_protein_cog_and_gene(ptt_file):
    """Build a dictionary cog_mapping PID to COG based on protein table file."""
    cog_mapping = {}
    product_mapping = {}
    with open(ptt_file) as read_handle:
        # Skip past either the first line, or first two lines, as NCBI can't seem to stick to any single one convention
        while True:
            # Line 1: Escherichia coli 536, complete genome - 1..4938920
            # Line 2: 4619 proteins
            line = read_handle.readline().strip()
            if line.endswith(' proteins'):
                break

        # Line 3: Location    Strand    Length    PID    Gene    Synonym    Code    COG    Product
        line = read_handle.readline().strip()
        columns = line.split("\t")
        assert columns == ['Location', 'Strand', 'Length', 'PID', 'Gene', 'Synonym', 'Code', 'COG', 'Product']

        for line in read_handle:
            values = line.strip().split('\t')
            assert len(values) == 9
            # Make format match format in genbank db_xref records
            pid = values[3]
            assert pid not in cog_mapping, 'Protein identifier was already assigned value: ' + cog_mapping[pid]
            assert pid not in product_mapping, 'Protein identifier was already assigned value: ' + product_mapping[pid]
            # Assign None if value of COG does not start with COG (such as when it's '-')
            cog_mapping[pid] = values[7] if values[7].startswith('COG') else None
            product_mapping[pid] = values[8].replace('<', '_').replace('>', '_').replace('|', '_').replace(';','_')
    return cog_mapping, product_mapping


def translate_genomes(genomes):
    """Download genome files, extract genes and translate those to proteins, returning DNA and protein fasta files."""
    assert len(genomes), 'Some genomes should be selected'

    # Use a pool to download files in the background while translating
    pool = Pool()
    futures = [pool.apply_async(download_genome_files, (genome,)) for genome in genomes]
    futures.extend([pool.apply_async(download_plasmid_files, (genome,)) for genome in genomes])
    dna_aa_pairs = [_translate_genome(gbk_ptt_pairs.get()) for gbk_ptt_pairs in futures if gbk_ptt_pairs.get() != None]

    # Extract DNA & Protein files separately from dna_aa_pairs
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

    # Concatenate files into one
    dna_concatemer = os.path.join(out_dir, '{pid}.ffn'.format(pid=project_id))
    protein_concatemer = os.path.join(out_dir, '{pid}.faa'.format(pid=project_id))
    concatenate(dna_concatemer, dna_files)
    concatenate(protein_concatemer, protein_files)
    return dna_concatemer, protein_concatemer


def _extract_gene_and_protein(out_dir, project_id, genbank_file, ptt_file=None, filetype=None):
    """Translate genbank DNA to protein and return resulting fasta files."""
    # Determine filenames for temporary and cache destination output files
    file_root = os.path.splitext(genbank_file)[0]
    dna_file_dest = os.path.join(out_dir, os.path.split(file_root)[1] + '.ffn')
    aa_file_dest = os.path.join(out_dir, os.path.split(file_root)[1] + '.faa')

    # Check if destination output files already exist and have content, and if so return them
    if (os.path.isfile(aa_file_dest) and 0 < os.path.getsize(aa_file_dest) and
        os.path.isfile(dna_file_dest) and 0 < os.path.getsize(dna_file_dest)):
        # But only if the output files are newer than the genbank file, otherwise translate the newer file & overwrite
        if (os.path.getmtime(genbank_file) < os.path.getmtime(aa_file_dest) and
            os.path.getmtime(genbank_file) < os.path.getmtime(dna_file_dest)):
            return dna_file_dest, aa_file_dest

    # Use temporary files as write handles when translating, so we can not pollute cache with incomplete files
    dna_tmp = tempfile.mkstemp(suffix='.ffn', prefix='translate_')[1]
    aa_tmp = tempfile.mkstemp(suffix='.faa', prefix='translate_')[1]

    # Determine if we'll add COG annotation from the ptt file
    if ptt_file is None:
        log.warn('No .ptt file given for %s: COG assignment disabled', genbank_file)
        # Disable COG mapping by setting cog_dict to an empty dictionary
        cog_dict = {}
        product_dict = {}
    else:
        # Retrieve PID to COG mapping from ptt file
        cog_dict, product_dict = _map_protein_cog_and_gene(ptt_file)

    # Open genbank_file & convert it using BioPython
    log.info('Translating %s', genbank_file)
    with open(aa_tmp, mode='w') as aa_wrtr:
        with open(dna_tmp, mode='w') as dna_wrtr:
            # Determine filetype by looking at extension, which should be either genbank or embl
            if filetype == None:
                filetype = os.path.splitext(genbank_file)[1][1:]
                if filetype == 'gbk':
                    filetype = 'genbank'

            # Bio.GenBank.Record
            gb_recrd = SeqIO.read(genbank_file, filetype)

            # Some RefSeq records such as 61583 do not contain nucleic acids, but rather refer to other files.
            # This means any extracted sequences will only consist of NN or XX, meaning we can't continue.
            str_seq = str(gb_recrd.seq)
            if re.match('^N+$', str_seq) or re.match('^X+$', str_seq):
                log.error('No nucleic acid sequence found in file %s, meaning we can not determine divergence for %s.',
                          genbank_file, project_id)

            # Determine whether the source of this genbank file is core genome or plasmid
            # plasmid = False
            # for gb_featr in gb_recrd.features:# Bio.SeqFeature
            #    if gb_featr.type == 'source' and 'plasmid' in gb_featr.qualifiers:
            #        plasmid = True
            #        break

            # Select only the coding sequences from all feature records
            coding_features = [gb_featr for gb_featr in gb_recrd.features
                               # Skip any non coding sequence features or pseudo (non-functional version) CDS
                               if (gb_featr.type == 'CDS'
                                   and not 'pseudo' in gb_featr.qualifiers
                                   and not 'pseudogene' in gb_featr.qualifiers
                                   )]

            # Remove duplicates while retaining order
            seen = set()
            coding_features = [x for x in coding_features if x not in seen and not seen.add(x)]

            # If there are no coding features, report this back to the user with a clear message rather than empty file
            if 0 == len(coding_features):
                log.error('No coding sequences found in file %s. Require protein table files to prevent this.',
                          genbank_file)

            # Translate all coding sequences
            for cds_featr in coding_features:
                _extract_and_translate_cds(cog_dict, product_dict, aa_wrtr, dna_wrtr, project_id, gb_recrd, cds_featr)

    assert os.path.isfile(dna_tmp) and 0 < os.path.getsize(dna_tmp), dna_tmp + ' should exist and have some content'
    assert os.path.isfile(aa_tmp) and 0 < os.path.getsize(aa_tmp), aa_tmp + ' should exist and have some content'

    # Move completed files to cache location only just now, so incomplete files do not pollute cache when raising errors
    shutil.move(dna_tmp, dna_file_dest)
    shutil.move(aa_tmp, aa_file_dest)

    return dna_file_dest, aa_file_dest


def _extract_and_translate_cds(cog_mapping, product_mapping, aa_writer, dna_writer, project_id, gb_record, gb_feature):
    """Extract DNA for Coding sequences, translate using GBK translation table and return dna & protein fasta files."""

    # Protein identifier is a property of the genbank feature
    protein_id = gb_feature.qualifiers['protein_id'][0]

    # Original sequence retrieved through BioPython 1.53+'s internal method
    extracted_seq = gb_feature.extract(gb_record.seq)  # Bio.Seq.Seq

    if len(extracted_seq) % 3:
        # Skip CDS feature if length of extracted_seq is not a multiple of three
        log.warn('Length of extracted coding sequence %s not a multiple of 3: %i\n%s',
                  protein_id, len(extracted_seq), extracted_seq)
        return

    # Some records do not contain nucleic acids, but rather refer to other files for their contigs: We can't handle that 
    str_seq = str(extracted_seq)
    if re.match('^X+$', str_seq) or re.match('^N+$', str_seq):
        # Skip CDS feature if length of extracted_seq is not a multiple of three
        log.warn('Extracted sequence only consists of X or N: %s', protein_id)
        return

    # Translation table is a property of the genbank feature
    transl_table = gb_feature.qualifiers['transl_table'][0]
    codon_table = CodonTable.unambiguous_dna_by_id.get(int(transl_table))

    # Set flag only when this CDS ends in a stop codon, so we can strip it off later, but do not strip non-stop-codons
    cds_has_stopcodon = str(extracted_seq[-3:]) in codon_table.stop_codons

    # Translate entire sequence as coding sequence using above translation table
    # Additional CodonTables are optionally available from Bio.Data.CodonTable
    try:
        # If transl_except is present, an alternative translation is offered in the GenBank record
        if 'transl_except' in gb_feature.qualifiers:
            # Fall back on GenBank translation whenever a transl_except record is found
            protein_seq = gb_feature.qualifiers['translation'][0]
        else:
            protein_seq = extracted_seq.translate(table=transl_table, cds=True)
    except TranslationError as err:
        # Try to recover from the following errors
        if ('First codon ' in str(err) and ' is not a start codon' in str(err)) or \
           ('Final codon ' in str(err) and ' is not a stop codon' in str(err)):
            # Retrieve GenBank provided reference translation
            provided_seq = gb_feature.qualifiers['translation'][0]
            # Calculate lengths for DNA & protein sequences
            dna_length = len(extracted_seq)
            aa_length = len(provided_seq)

            # Occasionally an incomplete amino end coding sequence is found, such as in: YP_004377721.1
            log.warning('Incomplete CDS found for %s from %s:\n%s', protein_id, gb_record.id, extracted_seq)

            # Calculate lengths for DNA & protein sequences
            if cds_has_stopcodon:
                # Correct check for the presence of a final stop codon
                dna_length -= 3

            # Assert lengths match, when stopcodon length is deducted
            if dna_length != 3 * aa_length:
                log.warn('Triple protein length %i did not match DNA length %i: record will be skipped',
                         3 * aa_length, dna_length)
                # Skip this CDS feature
                return

            protein_seq = provided_seq
            log.warning('Reverted to GenBank provided translation:\n%s', protein_seq)
        else:
            # Log some debug information before reraising error
            log.warn('Error in translating %s from %s using:\n%s', protein_id, gb_record.id, extracted_seq)
            log.warn(gb_feature)
            log.warn(str(err) + ': record will be skipped')
            # Skip this CDS feature, but do continue translating future CDS within the same GenBank record
            return

    # Determine COG by looking it up based on protein identifier
    cog = None
    product = None
    idlist = []
    if 'db_xref' in gb_feature.qualifiers:  # Prevent key error that might occur with some EMBL records
        idlist = [xref[3:] for xref in gb_feature.qualifiers['db_xref'] if xref.startswith('GI:')]
    if idlist:
        pid = idlist[0]
        cog = cog_mapping.get(pid, None)
        product = product_mapping.get(pid, None)

    # Strip stopcodon from DNA sequence here, so it's removed by the time we align & trim the sequences
    if cds_has_stopcodon:
        extracted_seq = extracted_seq[:-3]

    # Write out fasta. Header format as requested: >project_id|genbank_ac|protein_id|cog|gene name
    header = '{0}|{1}|{2}|{3}|{4}'.format(project_id, gb_record.id, protein_id, cog, product)
    _write_fasta(aa_writer, header, str(protein_seq))
    _write_fasta(dna_writer, header, str(extracted_seq))


def _write_fasta(write_handle, headerline, sequence):
    """Wrapper to write fasta header and sequence to file in correct format with new lines and limited line length."""
    # Write fasta to output file
    write_handle.write('>' + headerline + '\n')
    # Print rows with a length of 70 characters
    while 70 < len(sequence):
        write_handle.write(sequence[:70] + '\n')
        sequence = sequence[70:]
    # Write any remainder
    if sequence:
        write_handle.write(sequence + '\n')


def translate_fasta_coding_regions(nucl_fasta_file):
    """Translate an individual nucleotide fasta file containing coding regions to proteins using NCBI codon table 11."""
    # Determine output file name
    record_iter = SeqIO.parse(nucl_fasta_file, 'fasta')
    genomeid = record_iter.next().id.split('|')[0]  # pylint: disable=E1101
    prot_fasta_file = tempfile.mkstemp(suffix='.faa', prefix=genomeid + '.')[1]
    with open(prot_fasta_file, mode='w') as write_handle:
        for nucl_seqrecord in SeqIO.parse(nucl_fasta_file, 'fasta', alphabet=ambiguous_dna):
            # Translate nucl_seqrecord.seq
            try:
                prot_sequence = nucl_seqrecord.seq.translate(table=BACTERIAL_CODON_TABLE)
            except TranslationError as trer:
                log.warn('Skipping sequence because of translation error:\n%s', nucl_seqrecord)
                log.warn(trer)
                continue

            # Create protein sequence record and write it to file
            prot_seqrecord = SeqRecord(prot_sequence, id=nucl_seqrecord.id, description='')
            SeqIO.write(prot_seqrecord, write_handle, 'fasta')
    return prot_fasta_file


def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: translate.py
--genomes=FILE         file with selected genome IDs followed by Organism Name on each line
--external-zip=FILE    optional archive of user provided external genomes containing formatted nucleotide fasta files
--dna-zip=FILE         destination file path for zip archive of extracted DNA files
--protein-zip=FILE     destination file path for zip archive of translated protein files
"""
    options = ['genomes', 'external-zip=?', 'dna-zip', 'protein-zip']
    genome_ids_file, external_zip, dna_zipfile, protein_zipfile = parse_options(usage, options, args)

    dna_files = []
    protein_files = []

    # Read GenBank Project IDs from genomes_file, each on their own line
    with open(genome_ids_file) as read_handle:
        genome_ids = [line.split()[0] for line in read_handle
                      if not line.startswith('#') and 'external genome' not in line]

        if len(genome_ids):
            # Retrieve associated genome dictionaries from complete genomes table
            genomes = select_genomes_by_ids(genome_ids).values()
            genomes = sorted(genomes, key=itemgetter('Organism/Name'))

            # Actually translate the genomes to produced a set of files for both  dna files & protein files
            dna_files, protein_files = translate_genomes(genomes)

    # Also translate the external genomes
    if external_zip:
        # Extract external genomes archive
        external_dir = tempfile.mkdtemp(prefix='external_genomes_')
        external_dna_files = extract_archive_of_files(external_zip, external_dir)

        # Append IDs of external fasta files to genome IDs file
        _append_external_genomes(external_dna_files, genome_ids_file)

        # Translate individual files
        external_protein_files = [translate_fasta_coding_regions(dna_file) for dna_file in external_dna_files]

        # Add the files to the appropriate collections
        dna_files.extend(external_dna_files)
        protein_files.extend(external_protein_files)

    # Write the produced files to command line argument filenames
    create_archive_of_files(dna_zipfile, dna_files)
    create_archive_of_files(protein_zipfile, protein_files)

    # Do not clean up extracted DNA files or Protein translations: Keep them as cache

    # But do clean up external_dir now that the compressed archives are created
    if external_zip:
        shutil.rmtree(external_dir)

    # Exit after a comforting log message
    log.info("Produced: \n%s &\n%s", dna_zipfile, protein_zipfile)

if __name__ == '__main__':
    main(sys.argv[1:])
