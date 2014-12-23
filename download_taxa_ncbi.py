#!/usr/bin/env python
"""Module to download files from NCBI FTP."""
from Bio import SeqIO
from datetime import datetime, timedelta
from ftplib import FTP, error_perm
import logging
import os.path
import shutil
import tempfile
import time

from shared import create_directory


__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


def download_plasmid_files(genome):
    '''
    Download plasmid files from FTP site, by taking accession codes from other columns.
    :param genome:
    '''
    return download_genome_files(genome, download_log=None, refseq_column='Plasmids/RefSeq')


def download_genome_files(genome, download_log=None, require_ptt=False, refseq_column='Chromosomes/RefSeq'):
    """Download genome .gbk & .ptt files from ncbi ftp and return pairs per accessioncode in tuples of three."""
    logging.debug('Downloading: %s', genome)

    # ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria/Sulfolobus_islandicus_M_14_25_uid18871/CP001400.ffn
    # Download using FTP
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login(passwd='timtebeek+odose@gmail.com')

    # Try to find project directory in RefSeq curated listing
    projectid = genome['Assembly Accession']
    folder = '/genomes/ASSEMBLY_BACTERIA/{}'.format(genome['FTP Path'])
    accessioncodes = genome[refseq_column]
    target_dir = create_directory('genomes/' + projectid)

    # Determine last modified date to see if we should redownload the file following changes
    last_change_date = genome['Modify Date'] if genome['Modify Date'] else genome['Release Date']

    # Download .gbk & .ptt files for all genome accessioncodes and append them to this list as tuples of gbk + ptt
    genome_files = []

    for acc in accessioncodes:
        # Remove version suffixes to accessioncodes, such as NC_0012345.2
        acc = acc.split('.')[0]

        # Try genbank file, which is always required
        try:
            gbk_file = _download_genome_file(ftp, folder, acc + '.gbk', target_dir, last_change_date)

            # Try to parse Bio.GenBank.Record to see if it contains more than five (arbitrary) feature records
            features = SeqIO.read(gbk_file, 'genbank').features
            if not any(feature.type == 'CDS' for feature in features):
                # Skip when genbank file does not contain any coding sequence features
                logging.warn('GenBank file %s did not contain any coding sequence features', acc)
                continue
        except error_perm as err:
            if 'No such file or directory' not in str(err):
                raise err
            logging.warn(err)
            logging.warn('GenBank file %s missing for %s', acc, projectid)
            continue
        except IOError as err:
            if 'Target file was empty after download' not in str(err):
                raise err
            logging.warn(err)
            continue

        # Try protein table file, which could be optional
        ptt_file = None
        try:
            ptt_file = _download_genome_file(ftp, folder, acc + '.ptt', target_dir, last_change_date)
        except error_perm as err:
            if 'No such file or directory' not in str(err):
                raise err
            logging.warn(err)
            if require_ptt:
                logging.warn('Protein table file %s missing for %s: Probably no coding sequences', acc, projectid)
                continue
        except IOError as err:
            if 'Target file was empty after download' not in str(err):
                raise err
            logging.warn(err)
            continue
        genome_files.append((projectid, gbk_file, ptt_file))

    # Be nice and close the connection
    ftp.close()

    if len(genome_files) == 0:
        # Write out commented out line to the logfile detailing this error
        if download_log:
            with open(download_log, mode='a') as append_handle:
                append_handle.write('#{0}\t{1}\t'.format(projectid, genome['Organism/Name']))
                append_handle.write('# Genome skipped because of missing files\n')

        # Return nothing when:
        #- none of the accessioncodes resulted in files
        #- there were no protein table files when they were required
        #- no folder could be found for projectid
        return None

    # Write out provenance logfile with sources of retrieved files
    # This file could coincidentally also serve as genome ID file for extract taxa
    if download_log:
        with open(download_log, mode='a') as append_handle:
            append_handle.write('{0}\t{1}\t{2}{3}\n'.format(projectid, genome['Organism/Name'], ftp.host, folder))

    # Extend with
    if refseq_column == 'Chromosomes/RefSeq':
        genome_files.extend(download_plasmid_files(genome))

    # Return genome files
    return genome_files


def _find_project_dir(ftp, base_dir, projectid):
    """Find a genome project directory in ftp directory based upon directory name postfix. Return None if not found."""
    # Retrieve listing of directories under base_dir
    directory_listing = ftp.nlst(base_dir)

    # Find projectid dir based on directory suffix
    dir_suffix = '_uid{0}'.format(projectid)
    for ftp_file in directory_listing:
        if ftp_file.endswith(dir_suffix):
            return ftp_file
    return None


def _download_genome_file(ftp, remote_dir, filename, target_dir, last_change_date):
    """Download a single file from remote folder to target folder, only if it does not already exist."""
    # Move completed tmp_file to actual output path when done
    out_file = os.path.join(target_dir, filename)

    # We know when genomes were last updated. Use this information to determine when to download again, or every 60 days
    last_changed_stamp = time.mktime(last_change_date.timetuple())
    sixty_day_stamp = time.mktime((datetime.now() - timedelta(days=60)).timetuple())
    file_age_limit = max(last_changed_stamp, sixty_day_stamp)

    # Do not retrieve existing files if they exist and have content and are newer than the file_age_limit
    if not os.path.exists(out_file) or 0 == os.path.getsize(out_file) or os.path.getmtime(out_file) < file_age_limit:
        # Use a temporary file as write handle here, so we can not pollute cache when download is interrupted
        tmp_file = tempfile.mkstemp(prefix=filename + '_')[1]

        # Retrieve genbank & protein table files from FTP
        logging.info('Retrieving genome file %s%s/%s to %s', ftp.host, remote_dir, filename, target_dir)

        # Write retrieved contents to file
        with open(tmp_file, mode='wb') as writer:
            ftp.retrbinary('RETR {0}/{1}'.format(remote_dir, filename), writer.write)

        # Assert file was actually written to
        if not os.path.isfile(tmp_file) or 0 == os.path.getsize(tmp_file):
            raise IOError('Target file was empty after download: Did source have content?\n' + tmp_file)

        # Actually move now that we've finished downloading files
        shutil.move(tmp_file, out_file)
    else:
        logging.info('Cache hit on file %s dated %s', out_file, datetime.fromtimestamp(os.path.getmtime(out_file)))

    return out_file
