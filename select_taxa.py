#!/usr/bin/env python
"""Module for the select taxa step."""

from datetime import datetime, timedelta
from divergence import create_directory, HTTP_CACHE, parse_options, extract_archive_of_files
from ftplib import FTP, error_perm
from operator import itemgetter
import logging as log
import os
import shutil
import sys
import tempfile
import time

def select_genomes_by_ids(genome_ids):
    """Return list of genomes from complete genomes table whose GenBank Project ID is in genome_ids."""
    #Loop over genomes and return any genomes whose GenBank Project ID is in genome_ids
    refseq_genomes = dict((genome['RefSeq project ID'], genome) for genome in _parse_genomes_table())
    genbank_genomes = dict((genome['Project ID'], genome) for genome in _parse_genomes_table())

    #Require that the RefSeq and GenBank project IDs do not overlap, otherwise this fuzzy matching wont work anymore!
    assert set(refseq_genomes).isdisjoint(set(genbank_genomes)), 'RefSeq & GenBank project IDs should not overlap'

    #Match genomes_ids to genomes 
    matches = dict((queryid, refseq_genomes[queryid]) for queryid in genome_ids if queryid in refseq_genomes)
    matches.update(dict((queryid, genbank_genomes[queryid]) for queryid in genome_ids if queryid in genbank_genomes))

    #See if we matched all genomes, else log a warning
    for queryid in genome_ids:
        if queryid not in matches:
            log.warning('Could not find genome with Genbank Project ID %s in complete genomes table', queryid)

    return matches

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

def _parse_genomes_table(complete_genome_table = _download_genomes_table(), require_refseq = False):
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

        #Convert date columns to actual dates
        #Released date 01/27/2009
        if genome['Released date']:
            genome['Released date'] = datetime.strptime(genome['Released date'], '%m/%d/%Y')
        #Modified date 02/10/2011 17:38:40
        if genome['Modified date']:
            genome['Modified date'] = datetime.strptime(genome['Modified date'].split()[0], '%m/%d/%Y')

        #Append genome to list of genomes
        genomes.append(genome)

    #Filter out records not containing a refseq entry
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
    """Get tuples of Organism Name, GenBank Project ID & False, for input into Galaxy clade selection."""
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
                    name = '<b>{3}</b> - {0} &gt; {1} &gt; <i>{2}</i>'.format(
                        by_group(genome), by_firstname(genome), genome['Organism Name'], genome['Project ID'])
                    name += _get_colored_labels(genome)

                    #Yield the composed name, the project ID & False according to the expected input for Galaxy
                    yield name, genome['Project ID'], False

def _get_colored_labels(genome):
    """Optionally return colored labels for genome based on a release date, modified date and genome size."""
    labels = ''

    #Add New! & Updated! labels by looking at release and updated dates in the genome dictionary
    day_limit = datetime.today() - timedelta(days = 30)

    #Released date 01/27/2009
    released_date = genome['Released date']
    if released_date and day_limit < released_date:
        title = 'title="Since {0}"'.format(released_date.date())
        labels += ' <span {0} style="background-color: lightgreen">New!</span>'.format(title)

    #Modified date 02/10/2011 17:38:40
    modified_date = genome['Modified date']
    if modified_date and day_limit < modified_date:
        title = 'title="Since {0}"'.format(modified_date.date())
        labels += ' <span {0} style="background-color: yellow">Updated!</span>'.format(title)

    #Warn when genomes contain less than 0.5 MegaBase: unlikely to result in any orthologs
    genome_size = genome['Genome Size']
    if genome_size and float(genome_size) < 1:
        ttl = 'title="Small genomes are unlikely to result in orthologs present across all genomes"'
        style = 'style="background-color: orange"'
        labels += ' <span {0} {1}>Only {2} Mb!</span>'.format(ttl, style, genome_size)

    return labels

def download_genome_files(genome, download_log = None, require_ptt = False):
    """Download genome .gbk & .ptt files from ncbi ftp and return pairs per accessioncode in tuples of three."""
    #ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria/Sulfolobus_islandicus_M_14_25_uid18871/CP001400.ffn
    #Download using FTP
    ftp = FTP('ftp.ncbi.nih.gov')
    ftp.login()

    #Try to find project directory in RefSeq curated listing
    projectid = genome['RefSeq project ID']
    base_dir = '/genomes/Bacteria'
    project_dir = _find_project_dir(ftp, base_dir, projectid)
    if project_dir:
        accessioncodes = genome['List of RefSeq accessions']
        target_dir = create_directory('refseq/' + projectid)
    else:
        if projectid:
            log.warn('Genome directory not found under %s%s for %s', ftp.host, base_dir, projectid)

        #Try instead to find project directory in GenBank originals listing
        projectid = genome['Project ID']
        base_dir = '/genbank/genomes/Bacteria'
        project_dir = _find_project_dir(ftp, base_dir, projectid)
        if project_dir:
            accessioncodes = genome['List of GenBank accessions']
            target_dir = create_directory('genbank/' + projectid)
        else:
            log.warn('Genome directory not found under %s%s for %s', ftp.host, base_dir, projectid)

    #Determine ast modified date to see if we should redownload the file following changes
    last_change_date = genome['Modified date'] if genome['Modified date'] else genome['Released date']

    #Download .gbk & .ptt files for all genome accessioncodes and append them to this list as tuples of gbk + ptt
    genome_files = []

    #Occasionally we can not find a folder, meaning we will have to skip this genome as well
    if project_dir:
        for acc in accessioncodes:
            #Try genbank file, which is always required
            try:
                gbk_file = _download_genome_file(ftp, project_dir, acc + '.gbk', target_dir, last_change_date)
            except error_perm as err:
                if 'No such file or directory' not in str(err):
                    raise err
                log.warn(err)
                log.warn('GenBank file %s missing for %s', acc, projectid)
                continue

            #Try protein table file, which could be optional
            try:
                ptt_file = _download_genome_file(ftp, project_dir, acc + '.ptt', target_dir, last_change_date)
            except error_perm as err:
                if 'No such file or directory' not in str(err):
                    raise err
                log.warn(err)
                if require_ptt:
                    log.warn('Protein table file %s missing for %s: Probably no coding sequences', acc, projectid)
                    continue
                else:
                    ptt_file = None
            genome_files.append((projectid, gbk_file, ptt_file))

    #Be nice and close the connection
    ftp.close()

    if len(genome_files) == 0:
        #Write out commented out line to the logfile detailing this error
        if download_log:
            with open(download_log, mode = 'a') as append_handle:
                append_handle.write('#{0}\t{1}\t'.format(projectid, genome['Organism Name']))
                append_handle.write('#Genome skipped because of missing files\n')

        #Return nothing when:
        #- none of the accessioncodes resulted in files
        #- there were no protein table files when they were required
        #- no folder could be found for projectid
        return None

    #Write out provenance logfile with sources of retrieved files
    #This file could coincidentally also serve as genome ID file for extract taxa
    if download_log:
        with open(download_log, mode = 'a') as append_handle:
            append_handle.write('{0}\t{1}\t{2}{3}\n'.format(projectid, genome['Organism Name'], ftp.host, project_dir))

    #Return genome files
    return genome_files

def _find_project_dir(ftp, base_dir, projectid):
    """Find a genome project directory in ftp directory based upon directory name postfix. Return None if not found."""
    #Retrieve listing of directories under base_dir
    directory_listing = ftp.nlst(base_dir)

    #Find projectid dir based on directory suffix
    dir_suffix = '_uid{0}'.format(projectid)
    for ftp_file in directory_listing:
        if ftp_file.endswith(dir_suffix):
            return ftp_file
    return None

def _download_genome_file(ftp, remote_dir, filename, target_dir, last_change_date):
    """Download a single file from remote folder to target folder, only if it does not already exist."""
    #Move completed tmp_file to actual output path when done
    out_file = os.path.join(target_dir, filename)

    #We know when genomes were last updated. Use this information to determine when to download again, or every 60 days
    last_changed_stamp = time.mktime(last_change_date.timetuple())
    sixty_day_stamp = time.mktime((datetime.now() - timedelta(days = 60)).timetuple())
    file_age_limit = max(last_changed_stamp, sixty_day_stamp)

    #Do not retrieve existing files if they exist and have content and are newer than the file_age_limit
    if not os.path.exists(out_file) or 0 == os.path.getsize(out_file) or os.path.getmtime(out_file) < file_age_limit:
        #Use a temporary file as write handle here, so we can not pollute cache when download is interrupted
        tmp_file = tempfile.mkstemp(prefix = filename + '_')[1]

        #Retrieve genbank & protein table files from FTP
        log.info('Retrieving genome file %s from %s%s to %s', filename, ftp.host, remote_dir, target_dir)

        #Write retrieved contents to file
        with open(tmp_file, mode = 'wb') as write_file:
            #Write contents of file to shared experiment sra_lite
            def download_callback(block):
                """ftp.retrbinary requires a callback method"""
                write_file.write(block)
            ftp.retrbinary('RETR {0}/{1}'.format(remote_dir, filename), download_callback)

        #Assert file was actually written to
        assert os.path.isfile(tmp_file) and 0 < os.path.getsize(tmp_file), 'File should have content now: ' + tmp_file

        #Actually move now that we've finished downloading files
        shutil.move(tmp_file, out_file)

    return out_file

def _parse_external_genomes(external_zip, genomes_file):
    """Read out user provided labels and original filenames for uploaded genomes and append them to genome IDs file."""
    external_dir = tempfile.mkdtemp(prefix = 'external_genomes_')
    #External genomes are nucleotide fasta files uploaded by the user of which we have formatted the header
    external_fasta_files = extract_archive_of_files(external_zip, external_dir)

    for fasta_file in external_fasta_files:
        #Read user supplied species label & originating filename
        with open(fasta_file) as read_handle:
            #Read first line, which should be header of first record, starting with a skipped >
            header = read_handle.readline()[1:]

            #Header format as requested: >project_id|genbank_ac|protein_id|cog|source 
            label, origin = header.split('|')[0:2]

            with open(genomes_file, mode = 'a') as append_handle:
                #We'll use this 'external genome' source to skip externally derived files when downloading & translating
                append_handle.write('{0}\t{1}\texternal genome\n'.format(label, origin))

    #Remove temporarily extracted files
    shutil.rmtree(external_dir)

    return len(external_fasta_files)

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: select_taxa.py
--genomes=ID,...           optional comma-separated list of selected GenBank Project IDs from complete genomes table
--previous-file=FILE       optional previously or externally created GenBank Project IDs file whose genomes should be reselected
--external-zip=FILE        optional archive of user provided external genomes containing formatted nucleotide fasta files
--require-protein-table    require protein table files to be present for all downloaded genomes
--genomes-file=FILE        destination path for file with selected genome IDs followed by Organism Name on each line
"""
    options = ['genomes=?', 'previous-file=?', 'external-zip=?', 'require-protein-table?', 'genomes-file']
    genomes_line, previous_file, external_zip, require_ptt, genomes_file = parse_options(usage, options, args)

    #Handle externally uploaded genomes
    nr_of_external_genomes = 0
    if external_zip:
        nr_of_external_genomes = _parse_external_genomes(external_zip, genomes_file)

    #Genome IDs selected by the user that refer to GenBank or RefSeq entries
    genome_ids = []

    #Split ids on comma
    if genomes_line:
        genome_ids.extend(val for val in genomes_line.split(',') if val)

    #Allow for input of previous or externally created genomes-file to rerun an analysis
    if previous_file:
        #Read previous GenBank Project IDs from previous_file, each on their own line
        with open(previous_file) as read_handle:
            genome_ids.extend(line.split()[0] for line in read_handle
                              #But skip external genomes as their IDs will fail to download 
                              if 'external genome' not in line)

    #Assert each clade contains enough IDs
    minimum = 2
    maximum = 100
    if not minimum <= len(genome_ids) + nr_of_external_genomes <= maximum:
        log.error('Expected no less than {0} and no more than {1} selected genomes'.format(minimum, maximum))
        sys.exit(1)

    #Retrieve genome dictionaries to get to Organism Name
    genomes = select_genomes_by_ids(genome_ids).values()
    genomes = sorted(genomes, key = itemgetter('Organism Name'))

    #Write IDs to file, with organism name as second column to make the project ID files more self explanatory.
    for genome in genomes:
        #Download files here, but ignore returned files: These can be retrieved from cache during extraction/translation
        download_genome_files(genome, genomes_file, require_ptt = require_ptt)

    #Exit after a comforting log message
    log.info("Produced: \n%s", genomes_file)

if __name__ == '__main__':
    main(sys.argv[1:])
