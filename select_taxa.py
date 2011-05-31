#!/usr/bin/env python
"""Module for the select taxa step."""

from divergence import create_directory, HTTP_CACHE, parse_options
from ftplib import FTP, error_perm
from operator import itemgetter
import logging as log
import os
import shutil
import sys
import tempfile
import time

def select_genomes_by_ids(genome_ids):
    """Select genomes from complete genomes table if their RefSeq ID is in genome_file and return them as list."""
    #Loop over genomes and return any genomes whose RefSeq project ID is in genome_ids
    genomes = [genome for genome in _parse_genomes_table() if genome['RefSeq project ID'] in genome_ids]

    #See if we matched all genomes, else log a warning
    selected_ids = [genome['RefSeq project ID'] for genome in genomes]
    for genome_id in genome_ids:
        if genome_id not in selected_ids:
            log.warning('Could not find genome with RefSeq Project ID %s in complete genomes table', genome_id)

    return genomes

def _download_genomes_table():
    """Download complete genome table over HTTP with caching, decode as iso-8859-1 and return string."""
    genome_table = os.path.join(create_directory('.'), 'complete-genomes-table.tsv')

    time_between_downloads = 24 * 60 * 60
    if not os.path.isfile(genome_table) or os.path.getmtime(genome_table) < time.time() - time_between_downloads:
        #Download complete genome table over HTTP with caching and decode as iso-8859-1
        url = 'http://www.ncbi.nlm.nih.gov/genomes/lproks.cgi?dump=selected&view=1&p1=5:0'
        content = HTTP_CACHE.request(url)[1].decode('iso-8859-1')
        with open(genome_table, mode = 'w') as write_handle:
            write_handle.write(content.encode('utf-8'))
    else:
        #Read previously saved complete genomes table file
        with open(genome_table) as read_handle:
            content = read_handle.read().decode('utf-8')
    return content

def _parse_genomes_table(complete_genome_table = _download_genomes_table(), require_refseq = True):
    """Parse table of genomes and return list of dictionaries with values per genome."""
    #Empty lists to hold column names and genome dictionaries
    columns = []
    genomes = []

    #Separators are listed as suffixes to column headers: We'll split their values accordingly
    splitable_columns = {}

    #Loop over lines in complete genome table downloaded from url or read from cache
    for line in complete_genome_table.split('\n'):
        #Ignore empty lines
        if len(line) == 0:
            continue
        #Split values on tab
        values = line.split('\t')
        #Ignore comments
        if values[0].startswith(('#', '"#')):
            #But do retrieve column names from second header line
            if values[0] == '## Columns:':
                #Strip "" from every column name
                columns = [x[1:-1] for x in values[1:]]

                #Build up dictionary of split-able headers and their separators
                for i, column in enumerate(columns):
                    for suffix, separator in {' (comma separated)': ',', ' (pipe separated)': '|'}.iteritems():
                        if column.endswith(suffix):
                            columns[i] = column[:-len(suffix)]
                            splitable_columns[columns[i]] = separator
            continue

        #Assert we got the correct number of values for every column
        msg = u'Expected number of columns {0} and values {1} to match:\n{2}'.format(len(columns), len(values), line)
        assert len(columns) == len(values), msg

        #Build up a dictionary mapping column names to values for each genome, and append it to the list of genomes
        genome = dict(zip(columns, values))

        #Split values based on separator mapping for column name in splitable_columns
        for column, value in genome.iteritems():
            if column in splitable_columns:
                genome[column] = [] if len(value) in (0, 1) else value.split(splitable_columns[column])
        genomes.append(genome)

        #Filter out record not containing a refseq entry
        if require_refseq:
            genomes = [genome for genome in genomes if genome['RefSeq project ID']]

    #Return the genome dictionaries
    return tuple(genomes)

def _bin_using_keyfunctions(genomes, keyfunctions):
    """Bin genomes recursively according to keyfunctions, returning nested dictionaries mapping keys to collections."""

    def _bin_according_to_keyfunction(genomes, keyfunction = lambda x: x):
        """Group genomes into lists by unique values for keyfunction. Return dictionary mapping keys to lists."""
        bins = {}
        #Loop over genomes to map each to the right group
        for genome in genomes:
            #Determine key from keyfunction
            key = keyfunction(genome)
            #Retrieve either previous collection for key, or empty list, and assign back to bins[key]
            bins[key] = bins.get(key, [])
            #Add this genome to this colllection
            bins[key].append(genome)

        #Sort resulting lists by keyfunction as well, so we get a nice ordering of results in the end
        bins = dict((key, sorted(col, key = keyfunction)) for key, col in bins.iteritems())
        return bins

    #Get first keyfunction, and bin genomes by that keyfunction
    bins = _bin_according_to_keyfunction(genomes, keyfunctions[0])

    #If there are further keyfunctions available, recursively bin groups identified by current bin function
    furtherfunctions = keyfunctions[1:]
    if furtherfunctions:
        for key, subset in bins.iteritems():
            bins[key] = _bin_using_keyfunctions(subset, furtherfunctions)

    #Return dictionary of keys in genomes mapped to lists of genomes for that key, or further nested dictionaries
    return bins

def get_complete_genomes(genomes = _parse_genomes_table()):
    """Get tuples of Organism Name, RefSeq project ID & False, for input into Galaxy clade selection."""
    #Bin genomes using the following key functions iteratively
    by_kingdom = itemgetter('Super Kingdom')
    by_group = itemgetter('Group')
    by_firstname = lambda x:x['Organism Name'].split()[0]
    first_bin = _bin_using_keyfunctions(genomes, [by_kingdom, by_group, by_firstname])

    for key1 in sorted(first_bin.keys()):
        dict2 = first_bin[key1]
        for key2 in sorted(dict2.keys()):
            dict3 = dict2[key2]
            for key3 in sorted(dict3.keys()):
                list4 = dict3[key3]
                for genome in list4:
                    name = '<b>{3}</b> - {0} &gt; {1} &gt; <i>{2}</i>'.format(\
                        by_group(genome), by_firstname(genome), genome['Organism Name'], genome['RefSeq project ID'])
                    yield name, genome['RefSeq project ID'], False

def download_genome_files(genome):
    """Download genome .gbk & .ptt files from ncbi ftp and return pairs per accessioncode in tuples."""
    #ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria/Sulfolobus_islandicus_M_14_25_uid18871/CP001400.ffn
    host = 'ftp.ncbi.nih.gov'

    #Download using FTP
    ftp = FTP(host)
    ftp.login()

    #Try to find project directory in RefSeq curated listing
    projectid = genome['RefSeq project ID']
    base_dir = '/genomes/Bacteria'
    project_dir = _find_project_dir(ftp, base_dir, projectid)
    if project_dir:
        accessioncodes = genome['List of RefSeq accessions']
        target_dir = create_directory('refseq/' + projectid)
    else:
        log.warn('Genome directory not found under %s%s for %s', host, base_dir, projectid)

        #Try instead to find project directory in GenBank originals listing
        projectid = genome['Project ID']
        base_dir = '/genbank/genomes/Bacteria'
        project_dir = _find_project_dir(ftp, base_dir, projectid)
        if project_dir:
            accessioncodes = genome['List of GenBank accessions']
            target_dir = create_directory('genbank/' + projectid)
        else:
            log.error('Genome directory not found under %s%s for %s', host, base_dir, projectid)

    assert project_dir, 'Failed to find folder for genome {0}: {1}'\
        .format(genome['RefSeq project ID'], genome['Organism Name'])

    #Download .gbk & .ptt files for all genome refseq accessioncodes and append them to this list as tuples of gbk + ptt
    genome_files = []
    for acc in accessioncodes:
        gbk_file = _download_genome_file(ftp, project_dir, acc + '.gbk', target_dir)
        try:
            ptt_file = _download_genome_file(ftp, project_dir, acc + '.ptt', target_dir)
        except error_perm as err:
            if 'No such file or directory' in str(err):
                log.warn(err)
                ptt_file = None
            else:
                raise err
        genome_files.append((gbk_file, ptt_file))

    #Be nice and close the connection
    ftp.close()

    #Return genome files
    return genome_files

def _find_project_dir(ftp, base_dir, projectid):
    """Find a genome project directory in ftp directory based upon directory name postfix. Return None if not found."""
    #Retrieve listing of directories under refseq_dir
    directory_listing = ftp.nlst(base_dir)

    #Find refseq dir based on directory suffix
    dir_suffix = '_uid{0}'.format(projectid)
    for ftp_file in directory_listing:
        if ftp_file.endswith(dir_suffix):
            return ftp_file
    return None

def _download_genome_file(ftp, remote_dir, filename, target_dir):
    """Download a single file from remote folder to target folder, only if it does not already exist."""
    #Use a temporary file as write handle here, so we can not pollute cache when download is interrupted
    tmp_file = tempfile.mkstemp(prefix = filename)[1]

    #Do not retrieve existing files
    if not os.path.exists(tmp_file) or 0 == os.path.getsize(tmp_file):
        #Retrieve genbank & protein table files from FTP
        log.info('Retrieving genome file %s from %s%s to %s', filename, ftp.host, remote_dir, target_dir)

        #Write retrieved contents to file
        with open(tmp_file, mode = 'wb') as write_file:
            #Write contents of file to shared experiment sra_lite
            def download_callback(block):
                """ftp.retrbinary requires a callback method"""
                write_file.write(block)
            ftp.retrbinary('RETR {0}/{1}'.format(remote_dir, filename), download_callback)
    assert os.path.isfile(tmp_file) and 0 < os.path.getsize(tmp_file), 'File should exist and include some content now'

    #Move completed tmp_file to actual output path
    out_file = os.path.join(target_dir, filename)
    shutil.move(tmp_file, out_file)
    return out_file

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: select_taxa.py 
--genomes         comma-separated list of selected RefSeq project IDs from complete genomes table
--genomes-file    destination path for file with selected RefSeq project IDs followed by Organism Name on each line
"""
    options = ['genomes', 'genomes-file']
    genomes_line, genomes_file = parse_options(usage, options, args)

    #Split clade_a_ids & clade_b_ids each on comma
    genome_ids = [val for val in genomes_line.split(',') if val]

    #Assert each clade contains enough IDs
    minimum = 2
    maximum = 40
    if not minimum <= len(genome_ids) <= maximum:
        log.error('Expected no less than {0} and no more than {1} selected genomes'.format(minimum, maximum))
        sys.exit(1)

    #Retrieve genome dictionaries to get to Organism Name
    genomes = select_genomes_by_ids(genome_ids)

    #Write IDs to file, with organism name as second column, to make the RefSeq id files more self explanatory. 
    with open(genomes_file, mode = 'w') as write_handle:
        for genome in genomes:
            write_handle.write('{0}\t{1}\n'.format(genome['RefSeq project ID'], genome['Organism Name']))

    #Exit after a comforting log message
    log.info("Produced: \n%s", genomes_file)
    return genomes_file

if __name__ == '__main__':
    main(sys.argv[1:])
