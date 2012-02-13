#!/usr/bin/env python
"""Module to download files from MRS."""
from Bio import SeqIO
from datetime import datetime, timedelta
from divergence import create_directory
import logging
import os.path
import time
import urllib2

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


def download_genome_files(genome, download_log=None, require_ptt=False):
    """
    Download genome .gbk & .ptt files from MRS and return tuples containing project, genbank file & ptt file per
    accessioncode.

    @param genome: dictionary with genome values as parsed from lproks.cgi
    @param download_log: download log to append a line to for this genome
    @param require_ptt: boolean to indicate if individual accessioncodes should be skipped when ptt file is missing
    """
    #Try RefSeq accessions
    if genome['List of RefSeq accessions']:
        project = genome['RefSeq project ID']
        accessioncodes = genome['List of RefSeq accessions']
        databank = 'refseq'
        ptt_available = True
    else:
        project = genome['Project ID']
        accessioncodes = genome['List of GenBank accessions']
        databank = 'embl'
        ptt_available = False

        #MRS only has ptt files for Refseq: Fail early
        if require_ptt:
            return None

    #TODO Figure out how MRS handles updates, and whether we should fall back on genbank when refseq is missing

    #Determine output directory
    output_dir = create_directory(os.path.join(databank, project))

    #Determine last modified date to see if we should redownload the file following changes
    last_change_date = genome['Modified date'] if genome['Modified date'] else genome['Released date']

    #Download .gbk & .ptt files for all accessioncodes and append them to this list as tuples of project, gbk & ptt
    genome_files = []

    #Download all gbk & ptt files
    for ac in accessioncodes:
        #TODO This currently contains no error handling what so ever, while identifiers might in some cases be too new
        genbank_file = _download_file(output_dir, databank, ac, last_change_date)

        #Try to parse Bio.GenBank.Record to see if it contains more than five (arbitrary) feature records
        filetype = os.path.splitext(genbank_file)[1][1:]
        features = SeqIO.read(genbank_file, filetype).features
        if not any(feature.type == 'CDS' for feature in features):
            #Skip when genbank file does not contain any coding sequence features
            logging.warn('GenBank file %s did not contain any coding sequence features', ac)
            continue

        ptt_file = None if not ptt_available else _download_file(output_dir, 'ptt', ac, last_change_date)

        #Skip this accession when required ptt file is missing, but to allow for other accessions to pass
        if require_ptt and ptt_file == None:
            logging.warn('Protein table file %s missing for %s: Probably no coding sequences', ac, project)
            continue

        #Append tuples to genome_files
        genome_files.append((project, genbank_file, ptt_file))

    if len(genome_files) == 0:
        #Write out commented out line to the logfile detailing this error
        if download_log:
            with open(download_log, mode='a') as append_handle:
                append_handle.write('#{0}\t{1}\t'.format(project, genome['Organism Name']))
                append_handle.write('#Genome skipped because of missing files\n')

        #Return nothing when:
        #- none of the accessioncodes resulted in files
        #- there were no protein table files when they were required
        return None

    #Write out provenance logfile with sources of retrieved files
    #This file could coincidentally also serve as genome ID file for extract taxa
    if download_log:
        with open(download_log, mode='a') as append_handle:
            append_handle.write('{0}\t{1}\t{2}\n'.format(project, genome['Organism Name'],
                                                         'http://mrs.cmbi.ru.nl/mrs-5/info?db=' + databank))

    return genome_files


BASE_URL = 'http://mrs.cmbi.ru.nl/mrs-5/download?db={db}&id={id}'


def _download_file(output_dir, databank, ac, last_change_date):
    """Download a single databank file by accessioncode from MRS"""

    #Determine target output file path
    out_file = os.path.join(output_dir, ac + '.' + databank)

    #We know when genomes were last updated. Use this information to determine when to download again, or every 60 days
    last_changed_stamp = time.mktime(last_change_date.timetuple())
    sixty_day_stamp = time.mktime((datetime.now() - timedelta(days=60)).timetuple())
    file_age_limit = max(last_changed_stamp, sixty_day_stamp)

    #Do not retrieve existing files if they exist and have content and are newer than the file_age_limit
    if not os.path.exists(out_file) or 0 == os.path.getsize(out_file) or os.path.getmtime(out_file) < file_age_limit:
        url = BASE_URL.format(db=databank, id=ac)
        logging.info('Retrieving genome file %s to %s', url, out_file)

        #Download file from MRS
        response = urllib2.urlopen(url, timeout=60)
        content = response.read()

        #FIXME Temporarily stripping off content up to header line due to bug in MRS: Remove when remote bug is fixed
        if databank == 'refseq':
            from itertools import dropwhile
            content = ''.join(dropwhile(lambda x: not x.startswith('LOCUS') and not x.startswith('ID   '),
                                        content.splitlines(True)))

        #Save to local path
        with open(out_file, mode='w') as write_handle:
            write_handle.write(content)

        #Assert file was actually written to
        if not os.path.isfile(out_file) or 0 == os.path.getsize(out_file):
            raise IOError('Target file was empty after download: Did source have content?\n' + url)
    else:
        logging.info('Cache hit on file %s dated %s', out_file, datetime.fromtimestamp(os.path.getmtime(out_file)))

    return out_file


if __name__ == '__main__':
    moddate = datetime.strptime('02/13/2012', '%m/%d/%Y')
    print _download_file('/tmp', 'refseq', 'nc_011753', moddate)
    print _download_file('/tmp', 'ptt', 'nc_011753', moddate)
