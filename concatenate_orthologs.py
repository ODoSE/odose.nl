#!/usr/bin/env python
"""Module to fold individual ortholog sequences back into large concatemers according to their source genome."""

from Bio import SeqIO
from divergence import create_directory, parse_options, extract_archive_of_files, create_archive_of_files
import logging as log
import os.path
import shutil
import sys
import tempfile

def concatemer_per_genome(run_dir, trimmed_sicos):
    """Create a concatemer DNA file per genome containing all aligned & trimmed SICO genes."""
    concatemer_dir = create_directory('concatemers', inside_dir = run_dir)
    log.info('Creating concatemers from {0} SICOs'.format(len(trimmed_sicos)))

    #Open trimmed concatemer write handles
    concatemer_files = []
    write_handles = {}

    #Loop over trimmed sico files to append each sequence to the right concatemer
    for trimmed_sico in trimmed_sicos:
        for seqr in SeqIO.parse(trimmed_sico, 'fasta'):
            #Sample header line: >58191|NC_010067.1|YP_001569097.1|COG4948MR|core                
            refseq_id = seqr.id.split('|')[0]

            #Try to retrieve write handle from dictionary of cached write handles per genome
            write_handle = write_handles.get(refseq_id)

            #If not found, create & store write handle on demand
            if not write_handle:
                #Build up output file path for trimmed SICO genes concatemer per genome
                concatemer_file = os.path.join(concatemer_dir, refseq_id + '.trim.concat.fna')
                concatemer_files.append(concatemer_file)

                #Open write handle
                write_handle = open(concatemer_file, mode = 'w')
                write_handles[refseq_id] = write_handle

                #Write initial fasta header
                write_handle.write('> {0}|trimmed concatemer\n'.format(refseq_id))

            #Write / Append sequence for ortholog to genome concatemer
            write_handle.write('{0}\n'.format(str(seqr.seq)))

    #Close genomes trimmed concatemer write handles 
    for write_handle in write_handles.values():
        write_handle.close()

    log.info('Created {0} concatemers from {1} SICOs'.format(len(concatemer_files), len(trimmed_sicos)))

    return sorted(concatemer_files)

def create_super_concatemer(concatemer_files, destination_path):
    """Concatenate individual genome concatemers into a single super-concatemer for easy import into MEGA viewer."""
    with open(destination_path, mode = 'w') as write_handle:
        for concatemer in concatemer_files:
            seqr = SeqIO.read(concatemer, 'fasta')
            SeqIO.write(seqr, write_handle, 'fasta')

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: concatenate_orthologs.py
--orthologs-zip=FILE      archive of orthologous genes in FASTA format
--concatemer-zip=FILE     destination file path for archive of concatemers per genome
--concatemer-file=FILE    destination file path for super-concatemer of all genomes
"""
    options = ['orthologs-zip', 'concatemer-zip', 'concatemer-file']
    orthologs_zip, target_concat_zip, target_concat_file = parse_options(usage, options, args)

    #Run filtering in a temporary folder, to prevent interference from simultaneous runs
    run_dir = tempfile.mkdtemp(prefix = 'concatenate_orthologs_')

    #Extract files from zip archive
    temp_dir = create_directory('orthologs', inside_dir = run_dir)
    ortholog_files = extract_archive_of_files(orthologs_zip, temp_dir)

    #Concatenate trimmed_files per genome
    concatemer_files = concatemer_per_genome(run_dir, ortholog_files)

    #Create super concatemer
    create_super_concatemer(concatemer_files, target_concat_file)

    #Create concatemer files archive on command line specified path
    create_archive_of_files(target_concat_zip, concatemer_files)

    #Remove unused files to free disk space 
    shutil.rmtree(run_dir)

    #Exit after a comforting log message
    log.info('Produced: \n%s\n%s', target_concat_zip, target_concat_file)
    return target_concat_zip, target_concat_file

if __name__ == '__main__':
    main(sys.argv[1:])
