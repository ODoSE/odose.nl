#!/usr/bin/env python
"""Module for the reciprocal blast step."""

from divergence import create_directory, concatenate
from divergence.versions import MAKEBLASTDB, BLASTN, BLASTP
from multiprocessing import Pool
from subprocess import check_call, STDOUT
import logging as log
import os

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


def reciprocal_blast(run_dir, good_proteins_fasta, fasta_files):
    """Create blast database for good_proteins_fasta, blast all fasta_files against this database & return hits."""
    #Create blast database, retrieve path & name
    db_dir, db_name = _create_blast_database(run_dir, good_proteins_fasta)

    #Blast individual fasta files against the made blast databank, instead of the much larger good_proteins_fasta
    pool = Pool()  # Pool size is initialized to multiprocessing.cpu_count(), or 1 if that fails.

    #Submit job for each fasta files, and store future result handles
    submitted_jobs = [pool.apply_async(_blast_file_against_database, (db_dir, db_name, fasta)) for fasta in fasta_files]
    #Retrieve blast result files from future result handles and store in x_vs_all_hits
    x_vs_all_hits = (job.get() for job in submitted_jobs)

    #Concatenate the individual blast result files into one
    allvsall = os.path.join(db_dir, 'all-vs-all.tsv')
    concatenate(allvsall, x_vs_all_hits)
    return allvsall


def _create_blast_database(run_dir, fasta_file, nucleotide=False):
    """Create blast database"""
    assert os.path.exists(MAKEBLASTDB) and os.access(MAKEBLASTDB, os.X_OK), 'Could not find or run ' + MAKEBLASTDB

    dbtype = 'nucl' if nucleotide else 'prot'
    db_dir = create_directory('blast', inside_dir=run_dir)
    db_name = 'my_{0}_blast_db'.format(dbtype)
    log_file = os.path.join(db_dir, 'makeblastdb.log')
    with open(log_file, mode='w') as open_file:
        command = [MAKEBLASTDB,
                   '-in', fasta_file,
                   '-dbtype', dbtype,
                   '-out', os.path.join(db_dir, db_name)]
        log.info('Executing: %s', ' '.join(command))
        check_call(command, stdout=open_file)
    return db_dir, db_name


def _blast_file_against_database(db_dir, blast_db, fasta_file, nucleotide=False):
    """Blast all genes from genomes one and two against all genomes"""
    blast_program = BLASTN if nucleotide else BLASTP
    assert os.path.exists(blast_program) and os.access(blast_program, os.X_OK), 'Could not find or run ' + blast_program

    #Determine output file name
    basename = os.path.splitext(os.path.split(fasta_file)[1])[0]
    hits_file = os.path.join(db_dir, basename + '-vs-all.tsv')

    #Actually run the blast program
    command = [blast_program,
               '-db', blast_db,
               '-query', fasta_file,
               '-outfmt', str(6),
               '-out', hits_file]
    log.info('Executing: %s', ' '.join(command))
    check_call(command, cwd=db_dir, stdout=open('/dev/null', mode='w'), stderr=STDOUT)

    #Sanity check
    assert os.path.isfile(hits_file) and 0 < os.path.getsize(hits_file), hits_file + ' should exist with some content'
    return hits_file
