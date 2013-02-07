#!/usr/bin/env python
"""Module for the reciprocal blast step using the SARA Life Science Grid Portal."""

from datetime import timedelta
from divergence import concatenate
from divergence.lsgp_client import run_application, submit_application_run, retrieve_run_result, send_request
from divergence.versions import LSGP_MAKEBLASTDB, LSGP_BLASTN, LSGP_BLASTP
import os
import shutil
import tempfile

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


def reciprocal_blast(good_proteins_fasta, fasta_files):
    """Create blast database for good_proteins_fasta, blast all fasta_files against this database & return hits."""
    # Create blast database, retrieve path & name
    database_url = _create_blast_database(good_proteins_fasta)

    # Submit job for each fasta files, and store jobid & hits file name tuples
    jobids_and_hits_filenames = [_submit_blast_run(database_url, fasta_file) for fasta_file in fasta_files]

    # Retrieve all results, which should only take neglishably longer than waiting for the slowest results
    x_vs_all_hits = [_retrieve_blast_hits(jobid, hits_file) for jobid, hits_file in jobids_and_hits_filenames]

    # Concatenate the individual blast result files into one
    allvsall = tempfile.mkstemp(suffix='.tsv', prefix='all_vs_all')[1]
    concatenate(allvsall, x_vs_all_hits)

    # Clean up local and remote
    for x_vs_all in x_vs_all_hits:
        os.remove(x_vs_all)
    #Remote database_url port :444 requires certificate authentication, which we remove to use the basic authentication
    database_url = database_url.replace(':444', '')
    send_request(database_url, method='DELETE')

    return allvsall


def _create_blast_database(fasta_file, nucleotide=False):
    """Create blast database"""
    dbtype = 'nucl' if nucleotide else 'prot'
    db_name = '{0}_blast_db'.format(dbtype)

    #Run MAKEBLASTDB remotely
    params = {'dbtype': dbtype, 'out': db_name, 'portaldb': 1}
    files = {'in-file[]': fasta_file}
    db_dir = run_application(LSGP_MAKEBLASTDB,
                             params=params,
                             files=files,
                             max_duration=timedelta(days=1).total_seconds())

    #Upload database back to LSG Portal
    with open(os.path.join(db_dir, 'db_url.txt')) as read_handle:
        database_url = read_handle.read().strip()

    #Remove local database directory
    shutil.rmtree(db_dir)

    return database_url


def blast_file_against_database(database_url, fasta_file, nucleotide=False):
    """
    Blast file against database, and immediately retrieve and return the resulting BLAST hits file.
    @param database_url:
    @param fasta_file:
    @param nucleotide:
    """
    jobid, hits_file = _submit_blast_run(database_url, fasta_file, nucleotide)
    outside_path = _retrieve_blast_hits(jobid, hits_file)
    return outside_path


def _submit_blast_run(database_url, fasta_file, nucleotide=False):
    """
    Submit a BLAST run, returning the created jobid and filename for BLAST hits.
    @param database_url:
    @param fasta_file:
    @param nucleotide:
    """
    blast_app = LSGP_BLASTN if nucleotide else LSGP_BLASTP
    # Determine output file name
    hits_file = os.path.splitext(os.path.split(fasta_file)[1])[0] + '-vs-all.tsv'
    # Run BLAST[P/N] remotely
    params = {'db[]': database_url, 'out': hits_file, 'outfmt': 6}
    files = {'query': fasta_file}
    jobid = submit_application_run(blast_app,
                                   params=params,
                                   files=files)
    return jobid, hits_file


def _retrieve_blast_hits(jobid, hits_file):
    """
    Wait for job to complete, retrieve hits and clean up temporary files. Return the path to the created hits_file.
    @param jobid: jobid for which to retrieve the BLAST hits
    @param hits_file: name of the created BLAST hits file to retrieve
    """
    #Wait for job to complete and retrieve results
    results_dir = retrieve_run_result(jobid,
                                      max_duration=timedelta(hours=12).total_seconds())

    # Construct the likely path to the local file
    local_path = os.path.join(results_dir, hits_file)
    assert os.path.isfile(local_path) and 0 < os.path.getsize(local_path), local_path + ' should exist with content'

    # Move result file out of results_dir
    name, ext = os.path.splitext(os.path.split(local_path)[1])
    outside_path = tempfile.mkstemp(suffix=ext, prefix=name)[1]
    shutil.move(local_path, outside_path)

    # Remove results_dir, as we would lose directory reference otherwise
    shutil.rmtree(results_dir)

    return outside_path
