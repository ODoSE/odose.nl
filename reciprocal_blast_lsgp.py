#!/usr/bin/env python
"""Module for the reciprocal blast step using the SARA Life Science Grid Portal."""

from divergence.lsgp_client import run_application, upload_database
import os
import shutil

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"

MAKEBLASTDB = 'makeblastdb/2.2.25'
BLASTN = 'blastn/2.2.25'
BLASTP = 'blastp/2.2.25'


def _create_blast_database(fasta_file, nucleotide=False):
    """Create blast database"""
    dbtype = 'nucl' if nucleotide else 'prot'
    db_name = '{0}_blast_db'.format(dbtype)

    #Run MAKEBLASTDB remotely
    params = {'dbtype': dbtype,
              'out': db_name}
    files = {'in-file[]': fasta_file}
    db_dir = run_application(MAKEBLASTDB, params=params, files=files)

    #Upload database back to LSG Portal
    database = os.path.join(db_dir, db_name + '.tgz')
    database_url = upload_database(database)

    #Remove local database directory
    shutil.rmtree(db_dir)

    return database_url


def _blast_file_against_database(database_url, fasta_file, nucleotide=False):
    """Blast all genes from genomes one and two against all genomes"""
    blast_app = BLASTN if nucleotide else BLASTP

    #Determine output file name
    hits_file = os.path.splitext(os.path.split(fasta_file)[1])[0] + '-vs-all.tsv'

    #Run BLAST[P/N] remotely
    params = {'db[]': database_url,
              'out': hits_file,
              'outfmt': 6}
    files = {'query': fasta_file}
    db_dir = run_application(blast_app, params=params, files=files)

    #Construct the likely path to the local file
    local_path = os.path.join(db_dir, hits_file)
    assert os.path.isfile(local_path) and 0 < os.path.getsize(local_path), local_path + ' should exist with some content'

    return local_path

#TODO Convert
#def reciprocal_blast(run_dir, good_proteins_fasta, fasta_files):
#    """Create blast database for good_proteins_fasta, blast all fasta_files against this database & return hits."""
#    #Create blast database, retrieve path & name
#    db_dir, db_name = _create_blast_database(run_dir, good_proteins_fasta)
#
#    #Blast individual fasta files against the made blast databank, instead of the much larger good_proteins_fasta
#    pool = Pool()  # Pool size is initialized to multiprocessing.cpu_count(), or 1 if that fails.
#
#    #Submit job for each fasta files, and store future result handles
#    submitted_jobs = [pool.apply_async(_blast_file_against_database, (db_dir, db_name, fasta)) for fasta in fasta_files]
#    #Retrieve blast result files from future result handles and store in x_vs_all_hits
#    x_vs_all_hits = (job.get() for job in submitted_jobs)
#
#    #Concatenate the individual blast result files into one
#    allvsall = os.path.join(db_dir, 'all-vs-all.tsv')
#    concatenate(allvsall, x_vs_all_hits)
#    return allvsall
