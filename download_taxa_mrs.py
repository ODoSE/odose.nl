#!/usr/bin/env python
"""Module to download files from MRS."""
from Bio import SeqIO
from datetime import datetime, timedelta
from divergence import create_directory
from lxml import etree
import logging
import os.path
import time
import urllib2

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"

BASE_URL = 'http://mrs.cmbi.ru.nl/mrs-5/download?db={db}&id={id}'


def download_genome_files(genome, download_log=None, require_ptt=False):
    """
    Download genome .gbk & .ptt files from MRS and return tuples containing project, genbank file & ptt file per
    accessioncode.

    @param genome: dictionary with genome values as parsed from lproks.cgi
    @param download_log: download log to append a line to for this genome
    @param require_ptt: boolean to indicate if individual accessioncodes should be skipped when ptt file is missing
    """
    project = genome['BioProject ID']
    #Try RefSeq accessions
    if genome['Chromosomes/RefSeq']:
        accessioncodes = genome['Chromosomes/RefSeq']
        databank = 'refseq'
        ptt_available = True
    else:
        #Use embl accessions
        accessioncodes = genome['Chromosomes/INSDC']
        databank = 'embl'
        ptt_available = False

        #MRS only has ptt files for Refseq: Fail early
        if require_ptt:
            return None

    #The MRS embl & refseq databases receive updates daily, so there should be no need to fallback from one to the other

    #Determine output directory
    output_dir = create_directory(os.path.join(databank, project))

    #Determine last modified date to see if we should redownload the file following changes
    last_change_date = genome['Modify Date'] if genome['Modify Date'] else genome['Release Date']

    #Download .gbk & .ptt files for all accessioncodes and append them to this list as tuples of project, gbk & ptt
    genome_files = []

    #Download all gbk & ptt files
    for acc in accessioncodes:
        #Version numbers such as in "NC_009801.1" are sometimes appended to accession numbers: This is problematic
        acc = acc.split('.')[0]

        try:
            genbank_file = _download_file(output_dir, databank, acc, last_change_date)
        except IOError as ioerr:
            logging.warn('{0} file {1} missing for {2} because of: {3}'.format(databank, acc, project, str(ioerr)))
            continue

        #Try to parse Bio.GenBank.Record to see if it contains any CDS feature records
        filetype = os.path.splitext(genbank_file)[1][1:]
        features = SeqIO.read(genbank_file, filetype).features
        if not any(feature.type == 'CDS' for feature in features):
            #Skip when genbank file does not contain any coding sequence features
            logging.warn('GenBank file %s did not contain any coding sequence features', acc)
            continue

        if not ptt_available:
            ptt_file = None
        else:
            try:
                ptt_file = _download_file(output_dir, 'ptt', acc, last_change_date)
            except IOError as ioerr:
                logging.warn('{0} file {1} missing for {2} because of: {3}'.format(databank, acc, project, str(ioerr)))
                ptt_file = None

        #Skip this accession when required ptt file is missing, but do allow for other accessions to pass
        if require_ptt and ptt_file == None:
            logging.warn('Protein table file %s missing for %s: Probably no coding sequences', acc, project)
            continue

        #Append tuples to genome_files
        genome_files.append((project, genbank_file, ptt_file))

    if len(genome_files) == 0:
        #Write out commented out line to the logfile detailing this error
        if download_log:
            with open(download_log, mode='a') as append_handle:
                append_handle.write('#{0}\t{1}\t'.format(project, genome['Organism/Name']))
                append_handle.write('#Genome skipped because of missing files\n')

        #Return nothing when:
        #- none of the accessioncodes resulted in files
        #- there were no protein table files when they were required
        return None

    #Write out provenance logfile with sources of retrieved files
    #This file could coincidentally also serve as genome ID file for extract taxa
    if download_log:
        with open(download_log, mode='a') as append_handle:
            append_handle.write('{0}\t{1}\t{2}\n'.format(project, genome['Organism/Name'],
                                                         'http://mrs.cmbi.ru.nl/mrs-5/info?db=' + databank))

    return genome_files


def _download_file(output_dir, databank, acc, last_change_date):
    """Download a single databank file by accessioncode from MRS"""

    #Determine target output file path
    extension = 'genbank' if databank == 'refseq' else databank
    out_file = os.path.join(output_dir, acc + '.' + extension)

    #We know when genomes were last updated. Use this information to determine when to download again, or every 60 days
    #When no last updated date is given for a genome assume right now, to always download the file again
    last_changed_stamp = time.mktime(last_change_date.timetuple()) if last_change_date else time.time()
    sixty_day_stamp = time.mktime((datetime.now() - timedelta(days=60)).timetuple())
    file_age_limit = max(last_changed_stamp, sixty_day_stamp)

    #Do not retrieve existing files if they exist and have content and are newer than the file_age_limit
    if not os.path.exists(out_file) or 0 == os.path.getsize(out_file) or os.path.getmtime(out_file) < file_age_limit:
        url = BASE_URL.format(db=databank, id=acc)
        logging.info('Retrieving genome file %s to %s', url, out_file)

        #Download file from MRS
        response = urllib2.urlopen(url, timeout=60)
        content = response.read()

        #Assert file was actually written to
        if len(content) == 0:
            raise IOError('No content retrieved from download: Did source have content?\n' + url)

        #MRS does not raise 404 Not found for missing entries, but returns html content: Verify that didn't happen
        if 'mrs.cmbi.ru.nl' in content or '<html' in content:
            try:
                root = etree.fromstring(content, parser=etree.HTMLParser())
                reason = root.find('.//pre').text
                message = 'Error in downloading file; ' + reason
            except:
                message = 'Error in downloading file; Full response was:\n' + content
            raise IOError(message)

        #Only when above checks pass save to local path
        with open(out_file, mode='w') as write_handle:
            write_handle.write(content)
    else:
        logging.info('Cache hit on file %s dated %s', out_file, datetime.fromtimestamp(os.path.getmtime(out_file)))

    return out_file
