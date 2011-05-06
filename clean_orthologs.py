#!/usr/bin/env python
"""Module to clean up orthologs after OrthoMCL step."""

from __future__ import division
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from divergence import create_directory, extract_archive_of_files, create_archive_of_files
from divergence.select_taxa import select_genomes_from_file
from multiprocessing import Pool
from operator import itemgetter
from subprocess import check_call, STDOUT
import getopt
import logging as log
import os
import shutil
import sys
import tempfile

def trim_and_concat_sicos(genomes, dna_files, groups_file):
    """Invoke all cleanup operations sequentially and return the trimmed SICO files and their genome concatemers."""
    #Subdivide orthologs into groups
    shared_single_copy, shared_multi_copy, non_shared = _extract_shared_orthologs(genomes, groups_file)

    #Run cleanup in a temporary folder, to prevent interference from simultaneous runs
    run_dir = tempfile.mkdtemp(prefix = 'cleanup_run_')

    #Extract fasta files per orthologs
    sico_files, nr_of_seqs = _dna_file_per_sico(run_dir, dna_files, shared_single_copy)

    #Write statistics file
    stats_file = _write_statistics_file(run_dir, genomes, shared_single_copy, shared_multi_copy, non_shared, nr_of_seqs)

    #Look at assigned COGs among orthologs; filter those with multiple COGs & carry over COGs to sequences without COG  
    sico_files = _cog_based_filtering(sico_files, stats_file)

    #Align SICO DNA sequences at protein level.
    dna_alignments = _align_sicos(run_dir, sico_files)

    #Trim to retain equal length sequences with non-gap start & end codon
    trimmed_sico_files = _trim_alignments(run_dir, dna_alignments, stats_file)

    #Create concatemer of trimmed SICOs for each genome, resulting in equal length genome concatemers
    concatemers = _concatemer_per_genome(run_dir, genomes, trimmed_sico_files)

    #Create archives outside run_dir ahead of run_dir removal
    trimmed_zip = tempfile.mkstemp('.zip', 'trimmed_run_')[1]
    concatemer_zip = tempfile.mkstemp('.zip', 'concatemer_run_')[1]
    create_archive_of_files(trimmed_zip, trimmed_sico_files)
    create_archive_of_files(concatemer_zip, concatemers)

    #Move stats_file outside run_dir as well ahead of run_dir removal
    target_stats_file = tempfile.mkstemp('.txt', 'stats_run_')[1]
    shutil.move(stats_file, target_stats_file)

    #Remove run_dir to free disk space
    shutil.rmtree(run_dir)

    return trimmed_zip, concatemer_zip, target_stats_file

def _create_ortholog_dictionaries(groups_file):
    """Convert groups file into a list of ortholog dictionaries, which map refseq_id to their associated proteins."""
    #Sample line: group_5332: 58017|YP_219088.1 58191|YP_001572431.1 59431|YP_002149136.1
    ortholog_proteins_per_genome = []
    with open(groups_file) as read_handle:
        for line in read_handle:
            #Start at 1 to ignore incremental generated group_id
            remainder = line.split()[1:]
            proteins_per_genome = {}
            for ortholog in remainder:
                refseq_id, protein_id = ortholog.split('|')
                #Use dict().get(key, fallback_value) here to retrieve and assign valid array values for missing keys 
                proteins_per_genome[refseq_id] = proteins_per_genome.get(refseq_id, [])
                proteins_per_genome[refseq_id].append(protein_id)
            #Assign proteins per genome dictionary to orthologs per group as 
            ortholog_proteins_per_genome.append(proteins_per_genome)
    return ortholog_proteins_per_genome

def _extract_shared_orthologs(genomes, groups_file):
    """Filter orthologs to retain shared single and multiple copy orthologs from the collection of genomes."""
    log.info('Extracting shared orthologs for %d genomes from %s', len(genomes), groups_file)
    ortholog_proteins_per_genome = _create_ortholog_dictionaries(groups_file)

    #Retrieve RefSeq project IDs for the above genomes, which will be used to filter extracted orthologs later
    selected_genome_ids = []
    for genome in genomes:
        selected_genome_ids.append(genome['RefSeq project ID'])

    #Group orthologs into the following categories
    shared_multi_copy = []
    shared_single_copy = []
    non_shared_orthologs = []
    for prot_per_genomes in ortholog_proteins_per_genome:
        #Assume both shared and single copy
        is_shared_genome = True
        is_single_copy = True

        #Test validity of above boolean statements against proteins per genome 
        for selected_genome in selected_genome_ids:
            if selected_genome not in prot_per_genomes:
                #If any selected genome does not have proteins in this ortholog, it is not shared (no need to go on)
                is_shared_genome = False
                break
            elif 1 < len(prot_per_genomes[selected_genome]):
                #If any selected genome has more than one protein in this ortholog, it is not a single copy ortholog
                is_single_copy = False

        #Based on the now validated above boolean statements, optionally add proteins per genome to shared collections
        if is_shared_genome:
            if is_single_copy:
                shared_single_copy.append(prot_per_genomes)
            else:
                shared_multi_copy.append(prot_per_genomes)
        else:
            non_shared_orthologs.append(prot_per_genomes)

    #Assert all orthologs are binned
    sum_length = len(shared_single_copy) + len(shared_multi_copy) + len(non_shared_orthologs)
    assert len(ortholog_proteins_per_genome) == sum_length, 'All orthologs should fall into any one of three groups'

    #Return three collections of dictionaries mapping refseq_id to proteins for all orthologs
    return shared_single_copy, shared_multi_copy, non_shared_orthologs

def _write_statistics_file(run_dir, genomes, shared_single_copy, shared_multi_copy, partially_shared, nr_of_seqs):
    """Write out file with some basic statistics about the genomes, orthologs and size of shared core genome."""
    stats_file = os.path.join(run_dir, 'stats.txt')
    with open(stats_file, mode = 'w') as writer:
        #Some easy statistics about genomes and orthologs
        nr_genomes = len(genomes)
        nr_shared_sico = len(shared_single_copy)
        nr_shared_muco = len(shared_multi_copy)
        nr_part_shared = len(partially_shared)
        nr_orthologs = nr_shared_sico + nr_shared_muco + nr_part_shared

        #Determine number of ORFans
        nr_orfans = nr_of_seqs
        nr_orfans -= sum(len(genome_proteins) for genome_proteins in shared_single_copy)
        nr_orfans -= sum(len(genome_proteins) for genome_proteins in shared_multi_copy)
        nr_orfans -= sum(len(genome_proteins) for genome_proteins in partially_shared)

        #Write statistics to file        
        writer.write('{0:6}\tGenomes\n'.format(nr_genomes))
        writer.write('{0:6}\tGenes\n'.format(nr_of_seqs))
        writer.write('{0:6}\tOrphan genes\n\n'.format(nr_orfans))

        def perc(number):
            """Calculate a number as percentage of the number of orthologs"""
            return number / nr_orthologs
        writer.write('{0:6}\t({1:6.2%})\tShared single-copy orthologs\n'.format(nr_shared_sico, perc(nr_shared_sico)))
        writer.write('{0:6}\t({1:6.2%})\tShared multi-copy orthologs\n'.format(nr_shared_muco, perc(nr_shared_muco)))
        writer.write('{0:6}\t({1:6.2%})\tPartially shared orthologs\n'.format(nr_part_shared, perc(nr_part_shared)))
        writer.write('{0:6}\t({1:6.1%})\tTotal number of orthologs\n'.format(nr_orthologs, perc(nr_orthologs)))

    assert os.path.isfile(stats_file) and 0 < os.path.getsize(stats_file), stats_file + ' should exist with content.'
    return stats_file

def _dna_file_per_sico(run_dir, dna_files, shared_single_copy):
    """Create fasta files with all sequences per ortholog."""
    #Delete & create directory to remove any previously existing SICO files
    sico_dir = create_directory('sico', inside_dir = run_dir)

    #Loop over DNA files to extract SICO genes from each genome to file per SICO
    sico_files = set()
    number_of_sequences = 0
    for dna_file in dna_files:
        log.info('Extracting SICO genes from %s', dna_file)
        with open(dna_file) as read_handle:
            for record in SeqIO.parse(read_handle, 'fasta'):
                number_of_sequences += 1

                #Sample header line:    >58191|NC_010067.1|YP_001569097.1|COG4948MR|core
                #Corresponding ortholog: {'58191': ['YP_001569097.1'], ...}
                split_header = record.id.split('|')
                refseq_id = split_header[0]
                protein_id = split_header[2]

                #Loop over SICO's to write each gene from this DNA file to it's target SICO file
                for number, sico in enumerate(shared_single_copy):
                    #Create sico_file, with same filename formula as above
                    sico_file = os.path.join(sico_dir, 'sico_{0:06}.ffn'.format(number))
                    sico_files.add(sico_file)

                    #Write header & sequence to SICO file if protein ID matches SICO mapped protein ID for refseq_id
                    if protein_id == sico[refseq_id][0]:
                        #Append to the SICO files here to group the orthologs from various genomes in the same file
                        with open(sico_file, mode = 'a') as write_handle:
                            SeqIO.write([record], write_handle, 'fasta')
    return sorted(sico_files), number_of_sequences

#TODO Add option to filter out SICOs when any ortholog has been flagged as 'mobile element', 'phage' or 'IS element'

def _cog_based_filtering(sico_files, stats_file):
    """Inspect COGs for sequences marked as orthologs by OrthoMCL, and append some details about this to stats_file."""
    #Retrieve SICO to cog dictionaries of cog conflicts & transferable cog annotations, and list of SICOs missing cog
    cog_conflicts, cog_transferable, cog_missing = _group_cog_issues(sico_files)

    #Filter out orthologs containing more than one COG annotation
    sico_files = [sico for sico in sico_files if sico not in cog_conflicts.keys()]

    #Transfer COGs by overwriting sico_files with correct COG set
    for sico_file, cog in cog_transferable.iteritems():
        seqrecords = SeqIO.to_dict(SeqIO.parse(sico_file, 'fasta')).values()
        with open(sico_file, mode = 'w') as write_handle:
            for seqr in seqrecords:
                #Sample header line: >58191|NC_010067.1|YP_001569097.1|COG4948MR|core
                #Or for missing COG: >58191|NC_010067.1|YP_001569097.1|None|core
                split = seqr.id.split('|')
                assert split[3] in (cog, 'None'), 'COG should be either {0} or None, but was {1}'.format(cog, split[3])
                split[3] = cog
                seqr = SeqRecord(seqr.seq, id = '|'.join(split), description = '')
                SeqIO.write(seqr, write_handle, 'fasta')

    #Append statistics to stats file
    _append_cog_statistics(stats_file, cog_conflicts, cog_transferable, cog_missing)

    return sico_files

def _group_cog_issues(sico_files):
    """Find issues with COG assignments within SICO files by looking at COG conflicts, transferable and missing COGs."""
    cog_conflicts = {}
    cog_transferable = {}
    cog_missing = []
    for sico_file in sico_files:
        with open(sico_file) as read_handle:
            cogs = set()
            unassigned_cog_found = False
            for record in SeqIO.parse(read_handle, 'fasta'):
                #Sample header line: >58191|NC_010067.1|YP_001569097.1|COG4948MR|core
                cog = record.id.split('|')[3]
                if cog in cogs:
                    continue
                if cog == 'None':
                    unassigned_cog_found = True
                    continue
                cogs.add(cog)
            if 0 == len(cogs):
                cog_missing.append(sico_file)
            elif 1 == len(cogs):
                if unassigned_cog_found:
                    cog_transferable[sico_file] = cogs.pop()
            elif 1 < len(cogs):
                cog_conflicts[sico_file] = cogs
    return cog_conflicts, cog_transferable, cog_missing


def _append_cog_statistics(stats_file, cog_conflicts, cog_transferable, cog_missing):
    """Append COG statistics to stats_file"""
    with open(stats_file, mode = 'a') as append_handle:
        if cog_conflicts:
            msg = 'Multiple COGs found in {0} SICOs'.format(len(cog_conflicts))
            log.warn(msg)
            append_handle.write('\n' + msg + ':\n')
            for sico_file in sorted(cog_conflicts.keys()):
                cogs = cog_conflicts[sico_file]
                append_handle.write('{0}:\t{1}'.format(os.path.split(sico_file)[1], '\t'.join(cogs) + '\n'))
        if cog_transferable:
            msg = 'COGs transfered in {0} SICOs'.format(len(cog_transferable))
            log.info(msg)
            append_handle.write('\n' + msg + ':\n')
            for sico_file in sorted(cog_transferable.keys()):
                cog = cog_transferable[sico_file]
                append_handle.write('{0}:\t{1}'.format(os.path.split(sico_file)[1], cog + '\n'))
        if cog_missing:
            msg = 'No COGs found in {0} SICOs'.format(len(cog_missing))
            log.info(msg)
            append_handle.write('\n' + msg + ':\n')
            append_handle.write('\n'.join(os.path.split(sico_file)[1] for sico_file in cog_missing) + '\n')

def _align_sicos(run_dir, sico_files):
    """Align all SICO files given as argument in parallel and return the resulting alignment files."""
    log.info('Aligning {0} SICO genes using TranslatorX & muscle.'.format(len(sico_files)))
    #We'll multiplex this embarrassingly parallel task using a pool of workers
    tuples = [(run_dir, sico_file) for sico_file in sico_files]
    return Pool().map(_run_translatorx, tuples)

TRANSLATORX = '/projects/divergence/software/translatorx/translatorx_v1.1.pl'

def _run_translatorx((run_dir, sico_file), translation_table = '11'):
    """Run TranslatorX to create DNA level alignment file of protein level aligned DNA sequences within sico_file."""
    assert os.path.exists(TRANSLATORX) and os.access(TRANSLATORX, os.X_OK), 'Could not find or run ' + TRANSLATORX

    #Determine output file name
    sico_base = os.path.splitext(os.path.split(sico_file)[1])[0]
    alignment_dir = create_directory('alignments/' + sico_base, inside_dir = run_dir)

    #Target output file
    file_base = os.path.join(alignment_dir, sico_base)
    dna_alignment = file_base + '.nt_ali.fasta'

    #Actually run the TranslatorX program
    command = [TRANSLATORX,
               '-i', sico_file,
               '-c', translation_table,
               '-o', file_base]
    check_call(command, stdout = open('/dev/null', 'w'), stderr = STDOUT)

    msg = 'Alignment file should exist and have some content now: {0}'.format(dna_alignment)
    assert os.path.isfile(dna_alignment) and 0 < os.path.getsize(dna_alignment), msg
    return dna_alignment

def _trim_alignments(run_dir, dna_alignments, stats_file):
    """Trim all DNA alignments using _trim_alignment (singular), and calculate some statistics about the trimming."""
    log.info('Trimming {0} DNA alignments from first non-gap codon to last non-gap codon'.format(len(dna_alignments)))

    #Create directory here, to prevent race-condition when folder does not exist, but is then created by another process
    trimmed_dir = create_directory('trimmed', inside_dir = run_dir)

    #algn_perct_tpls = [_trim_alignment(trimmed_dir, ali) for ali in dna_alignments]
    tuples = [(trimmed_dir, dna_alignment) for dna_alignment in dna_alignments]
    algn_perct_tpls = Pool().map(_trim_alignment, tuples)
    algn_perct_tpls = sorted(algn_perct_tpls, key = itemgetter(1))

    remaining_percts = [tpl[1] for tpl in algn_perct_tpls]
    with open(stats_file, mode = 'a') as append_handle:
        append_handle.write('\n{0:6} sequence alignments trimmed\n'.format(len(algn_perct_tpls)))
        average_retained = sum(remaining_percts) / len(remaining_percts)
        append_handle.write('{0:6.2%} sequence retained on average\n'.format(average_retained))
        append_handle.write('10 least percentages retained from sequence alignments:\n')
        least = '\n'.join('{0}: {1:6.2%}'.format(os.path.split(align)[1], perc) for align, perc in algn_perct_tpls[:10])
        append_handle.write(least + '\n')

    trimmed_alignments = [tpl[0] for tpl in algn_perct_tpls]
    return sorted(trimmed_alignments)

def _trim_alignment((trimmed_dir, dna_alignment)):
    """Trim alignment to retain first & last non-gapped codons across alignment, and everything in between (+gaps!)."""
    with open(dna_alignment) as read_handle:
        #Read single alignment from fasta file
        alignment = AlignIO.read(read_handle, 'fasta')
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

        return trimmed_file, trimmed_length / alignment_length

def _concatemer_per_genome(run_dir, genomes, trimmed_sicos):
    """Create a concatemer DNA file per genome containing all aligned & trimmed SICO genes."""
    concatemer_dir = create_directory('concatemers', inside_dir = run_dir)
    log.info('Creating {0} concatemers from {1} SICOs'.format(len(genomes), len(trimmed_sicos)))

    #Open trimmed concatemer write handles
    concatemer_sico_files = []
    write_handles = {}

    #For each genome, open a file for the trimmed SICO genes concatemer
    for genome in genomes:
        refseq_id = genome['RefSeq project ID']

        #Build up output file path
        concatemer_file = os.path.join(concatemer_dir, refseq_id + '.trimmed.concatemer.fasta')
        concatemer_sico_files.append(concatemer_file)

        #Open write handle
        write_handle = open(concatemer_file, mode = 'w')
        write_handles[refseq_id] = write_handle

        #Write initial fasta header
        write_handle.write('> {0}|{1}|trimmed concatemer\n'.format(refseq_id, genome['Organism Name']))

    #Loop over trimmed sico files to append each sequence to the right concatemer
    for trimmed_sico in trimmed_sicos:
        with open(trimmed_sico) as read_handle:
            for seqr in SeqIO.parse(read_handle, 'fasta'):
                #Sample header line: >58191|NC_010067.1|YP_001569097.1|COG4948MR|core                
                refseq_id = seqr.id.split('|')[0]
                write_handles[refseq_id].write('{0}\n'.format(str(seqr.seq)))

    #Close genomes trimmed concatemer write handles 
    for write_handle in write_handles.values():
        write_handle.close()

    return concatemer_sico_files

def main(args):
    """Main function called when run from command line or as part of pipeline."""

    def _parse_options(args):
        """Use getopt to parse command line argument options"""

        def _usage():
            """Print _usage information"""
            print """
Usage: clean_orthologs.py 
--genomes=FILE           file with refseq id from complete genomes table on each line 
--dna-zip=FILE           zip archive of extracted DNA files
--groups=FILE            file listing groups of orthologous proteins
--trimmed-zip=FILE       destination file path for archive of aligned & trimmed single copy orthologous (SICO) genes
--concatemer-zip=FILE    destination file path for archive of SICO concatemer per genome
--stats=FILE             destination file path for SICO cleanup statistics file
"""

        options = ['genomes', 'dna-zip', 'groups', 'trimmed-zip', 'concatemer-zip', 'stats']
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

    genome_ids_file, dna_zip, groups_file, target_trimmed, target_concatemer, target_stats_path = _parse_options(args)

    #Parse file containing RefSeq project IDs & retrieve associated genome dictionaries from complete genomes table
    genomes = select_genomes_from_file(genome_ids_file)

    #Extract files from zip archive
    temp_dir = tempfile.mkdtemp()
    dna_files = extract_archive_of_files(dna_zip, temp_dir)

    #Actually run cleanup
    trimmed_zip, concatemer_zip, stats_file = trim_and_concat_sicos(genomes, dna_files, groups_file)

    #Move produced files to command line specified output paths
    shutil.move(trimmed_zip, target_trimmed)
    shutil.move(concatemer_zip, target_concatemer)
    shutil.move(stats_file, target_stats_path)

    #Remove unused files to free disk space 
    shutil.rmtree(temp_dir)

    #Exit after a comforting log message
    log.info("Produced: \n%s\n%s\n%s", trimmed_zip, concatemer_zip, target_stats_path)
    return trimmed_zip, concatemer_zip, target_stats_path

if __name__ == '__main__':
    main(sys.argv[1:])
