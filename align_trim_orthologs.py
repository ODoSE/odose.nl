#!/usr/bin/env python
"""Module to align and trim orthologs after the OrthoMCL step."""

from __future__ import division
from Bio import AlignIO
from divergence import create_directory, extract_archive_of_files, create_archive_of_files, parse_options
from divergence.versions import TRANSLATORX
from multiprocessing import Pool
from operator import itemgetter
from subprocess import check_call, STDOUT
import logging as log
import os
import shutil
import sys
import tempfile

def _align_sicos(run_dir, sico_files):
    """Align all SICO files given as argument in parallel and return the resulting alignment files."""
    log.info('Aligning {0} SICO genes using TranslatorX & muscle.'.format(len(sico_files)))
    #We'll multiplex this embarrassingly parallel task using a pool of workers
    tuples = [(run_dir, sico_file) for sico_file in sico_files]
    return Pool().map(_run_translatorx, tuples)

def _run_translatorx((run_dir, sico_file), translation_table = '11'):
    """Run TranslatorX to create DNA level alignment file of protein level aligned DNA sequences within sico_file."""
    assert os.path.exists(TRANSLATORX) and os.access(TRANSLATORX, os.X_OK), 'Could not find or run ' + TRANSLATORX

    #Determine output file name
    sico_base = os.path.splitext(os.path.split(sico_file)[1])[0]
    alignment_dir = create_directory('alignments/' + sico_base, inside_dir = run_dir)

    #Created output file
    file_base = os.path.join(alignment_dir, sico_base)
    dna_alignment = file_base + '.nt_ali.fasta'

    #Actually run the TranslatorX program
    command = [TRANSLATORX,
               '-i', sico_file,
               '-c', translation_table,
               '-o', file_base]
    check_call(command, stdout = open('/dev/null', 'w'), stderr = STDOUT)

    assert os.path.isfile(dna_alignment) and 0 < os.path.getsize(dna_alignment), \
        'Alignment file should exist and have some content now: {0}'.format(dna_alignment)
    return dna_alignment

def _trim_alignments(run_dir, dna_alignments, retained_threshold, stats_file):
    """Trim all DNA alignments using _trim_alignment (singular), and calculate some statistics about the trimming."""
    log.info('Trimming {0} DNA alignments from first non-gap codon to last non-gap codon'.format(len(dna_alignments)))

    #Create directory here, to prevent race-condition when folder does not exist, but is then created by another process
    trimmed_dir = create_directory('trimmed', inside_dir = run_dir)

    #Use Pool().map again to scale trimming out over multiple cores. This requires tuple'd arguments however
    trim_tpls = Pool().map(_trim_alignment, ((trimmed_dir, dna_alignment) for dna_alignment in dna_alignments))

    remaining_percts = [tpl[3] for tpl in trim_tpls]
    trimmed_alignments = [tpl[0] for tpl in trim_tpls if retained_threshold <= tpl[3]]
    misaligned = [tpl[0] for tpl in trim_tpls if retained_threshold > tpl[3]]

    #Write trim statistics to file in such a way that they're easily converted to a graph in Galaxy
    with open(stats_file, mode = 'w') as append_handle:
        msg = '{0:6} sequence alignments trimmed'.format(len(trim_tpls))
        log.info(msg)
        append_handle.write('#' + msg + '\n')

        average_retained = sum(remaining_percts) / len(remaining_percts)
        msg = '{0:5.2}% sequence retained on average overall'.format(average_retained)
        log.info(msg)
        append_handle.write('#' + msg + '\n')

        filtered = len(misaligned)
        msg = '{0:6} orthologs filtered as they retained less than {1}%'.format(filtered, str(retained_threshold))
        log.info(msg)
        append_handle.write('#' + msg + '\n')

        append_handle.write('# Trimmed file\tOriginal length\tTrimmed length\tPercentage retained\n')
        for tpl in sorted(trim_tpls, key = itemgetter(3)):
            append_handle.write(os.path.split(tpl[0])[1] + '\t')
            append_handle.write(str(tpl[1]) + '\t')
            append_handle.write(str(tpl[2]) + '\t')
            append_handle.write('{0:.2f}\n'.format(tpl[3]))

    return sorted(trimmed_alignments), sorted(misaligned)

def _trim_alignment((trimmed_dir, dna_alignment)):
    """Trim alignment to retain first & last non-gapped codons across alignment, and everything in between (+gaps!).
    
    Return trimmed file, original length, trimmed length and percentage retained as tuple"""
    #Read single alignment from fasta file
    alignment = AlignIO.read(dna_alignment, 'fasta')
    #print '\n'.join([str(seqr.seq) for seqr in alignment])

    #Total alignment should be just as long as first sequence of alignment
    alignment_length = len (alignment[0])

    #After using protein alignment only for CDS, all alignment lengths should be multiples of three 
    assert alignment_length % 3 == 0, 'Length not a multiple of three: {} \n{2}'.format(alignment_length, alignment)

    #Assert all codons are either full length codons or gaps, but not a mix of gaps and letters such as AA- or A--
    for index in range(0, alignment_length, 3):
        for ali in alignment:
            codon = ali.seq[index:index + 3]
            assert not ('-' in codon and str(codon) != '---'), '{0} at {1} in \n{2}'.format(codon, index, alignment)

    #Loop over alignment, taking 3 DNA characters each time, representing a single codon
    first_full_codon_start = None
    last_full_codon_end = None
    for index in range(0, alignment_length, 3):
        codon_concatemer = ''.join([str(seqr.seq) for seqr in alignment[:, index:index + 3]])
        if '-' in codon_concatemer:
            continue
        if first_full_codon_start is None:
            first_full_codon_start = index
        else:
            last_full_codon_end = index + 3

    #Create sub alignment consisting of all trimmed sequences from full alignment
    trimmed = alignment[:, first_full_codon_start:last_full_codon_end]
    trimmed_length = len(trimmed[0])
    assert trimmed_length % 3 == 0, 'Length not a multiple of three: {} \n{2}'.format(trimmed_length, trimmed)

    #Write out trimmed alignment file
    trimmed_file = os.path.join(trimmed_dir, os.path.split(dna_alignment)[1])
    with open(trimmed_file, mode = 'w') as write_handle:
        AlignIO.write(trimmed, write_handle, 'fasta')

    #Assert file now exists with content
    assert os.path.isfile(trimmed_file) and os.path.getsize(trimmed_file), \
        'Expected trimmed alignment file to exist with some content now: {0}'.format(trimmed_file)

    return trimmed_file, alignment_length, trimmed_length, trimmed_length / alignment_length * 100

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: filter_orthologs.py
--orthologs-zip=FILE           archive of orthologous genes in FASTA format
--retained-threshold=PERC      filter orthologs that retain less than PERC % of sequence after trimming alignment 
--aligned-zip=FILE             destination file path for archive of aligned orthologous genes
--misaligned-zip=FILE          destination file path for archive of misaligned orthologous genes per retained-threshold
--trimmed-zip=FILE             destination file path for archive of aligned & trimmed orthologous genes
--stats=FILE                   destination file path for ortholog trimming statistics file
"""
    options = ['orthologs-zip', 'retained-threshold', 'aligned-zip', 'misaligned-zip', 'trimmed-zip', 'stats']
    orthologs_zip, retained_threshold, aligned_zip, misaligned_zip, trimmed_zip, target_stats_path = \
        parse_options(usage, options, args)

    #Convert retained threshold to integer, so we can fail fast if argument value format was wrong
    retained_threshold = int(retained_threshold)

    #Run filtering in a temporary folder, to prevent interference from simultaneous runs
    run_dir = tempfile.mkdtemp(prefix = 'align_trim_')

    #Extract files from zip archive
    temp_dir = create_directory('orthologs', inside_dir = run_dir)
    sico_files = extract_archive_of_files(orthologs_zip, temp_dir)

    #Align SICOs so all sequences become equal length sequences
    trim_stats_file = os.path.join(run_dir, 'trim-stats.txt')
    aligned_files = _align_sicos(run_dir, sico_files)

    #Filter orthologs that retain less than PERC % of sequence after trimming alignment    
    trimmed_files, misaligned_files = _trim_alignments(run_dir, aligned_files, retained_threshold, trim_stats_file)

    #Create archives of files on command line specified output paths & move trim_stats_file
    create_archive_of_files(aligned_zip, aligned_files)
    create_archive_of_files(misaligned_zip, misaligned_files)
    create_archive_of_files(trimmed_zip, trimmed_files)
    shutil.move(trim_stats_file, target_stats_path)

    #Remove unused files to free disk space 
    shutil.rmtree(run_dir)

    #Exit after a comforting log message
    log.info('Produced: \n%s\n%s\n%s\n%s', aligned_zip, misaligned_zip, trimmed_zip, target_stats_path)

if __name__ == '__main__':
    main(sys.argv[1:])
