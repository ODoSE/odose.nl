#!/usr/bin/env python
"""Package divergence"""

from pkg_resources import resource_filename  # @UnresolvedImport  # pylint: disable=E0611
from zipfile import ZipFile, ZIP_DEFLATED, is_zipfile
import Bio
import getopt
import httplib2
import logging
import os
import shutil
import sys

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"

#Configure logging LOG_FORMAT
LOG_FORMAT = '%(levelname)s\t%(asctime)s %(module)s.%(funcName)s:%(lineno)d\t%(message)s'
LOG_DATE_FORMAT = '%H:%M:%S'

#Logs WARNING messages and anything above to sys.stdout
logging.basicConfig(level=logging.INFO, stream=sys.stdout, format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT)

#Log ERROR messages to stderr separately; these will fail a tool run in Galaxy
STDERR_HANDLER = logging.StreamHandler(sys.stderr)
STDERR_HANDLER.setLevel(logging.ERROR)
STDERR_HANDLER.setFormatter(logging.Formatter(fmt=LOG_FORMAT, datefmt=LOG_DATE_FORMAT))
logging.root.addHandler(STDERR_HANDLER)

#Require at least version 1.53 op BioPython
assert 1.54 <= float(Bio.__version__), 'BioPython version 1.54 or higher is required'

#Using the standard NCBI Bacterial, Archaeal and Plant Plastid Code translation table (11)
CODON_TABLE_ID = 11

#Base output dir
BASE_OUTPUT_PATH = '../divergence-cache/'


def create_directory(dirname, inside_dir=BASE_OUTPUT_PATH):
    """Create a directory in the default output directory, and return the full path to the directory.

    Return directory if directory already exists, raise error if file by that name already exists."""
    filename = os.path.join(inside_dir, dirname)
    #For non-absolute paths, get filename relative to this module
    if filename[0] != '/':
        filename = resource_filename(__name__, filename)

    #If file exists and is a directory, return the existing directory unaltered
    if os.path.exists(filename):
        if os.path.isdir(filename):
            return filename
        else:
            raise IOError('Could not create directory {0}\nA file with that name already exists.')
    else:
        os.makedirs(filename)
        return filename

#Initialize shared cache for files downloaded through httplib2
HTTP_CACHE = httplib2.Http(create_directory('.cache'))


def concatenate(target_path, source_files):
    """Concatenate arbitrary number of files into target_path by reading and writing in binary mode.

    WARNING: The binary mode implies new line characters will NOT be added in between files!"""
    with open(target_path, mode='wb') as write_handle:
        for source_file in source_files:
            shutil.copyfileobj(open(source_file, mode='rb'), write_handle)
    assert os.path.isfile(target_path) and 0 < os.path.getsize(target_path), target_path + ' should exist with content'


def create_archive_of_files(archive_file, file_iterable):
    """Write files in file_iterable to archive_file, using only filename for target path within archive_file."""
    zipfile_handle = ZipFile(archive_file, mode='w', compression=ZIP_DEFLATED)
    if len(file_iterable):
        for some_file in file_iterable:
            zipfile_handle.write(some_file, os.path.split(some_file)[1])
    else:
        logging.warn('No files in file_iterable: %s will be empty!', archive_file)
        zipfile_handle.writestr('empty', '')
    zipfile_handle.close()
    assert is_zipfile(archive_file), 'File should now have been a valid zipfile: ' + archive_file


def extract_archive_of_files(archive_file, target_dir):
    """Extract all files from archive to target directory, and return list of files extracted."""
    extracted_files = []
    read_handle = ZipFile(archive_file, mode='r')
    for zipinfo in read_handle.infolist():
        extracted_path = read_handle.extract(zipinfo, path=target_dir)
        extracted_files.append(extracted_path)
    read_handle.close()
    assert extracted_files, 'At least some files should have been read from ' + archive_file
    return extracted_files


def parse_options(usage, options, args):
    """Parse command line arguments in args. Options require argument by default; flags are indicated with '?' postfix.

    Parameters:
    usage -- Usage string detailing command line arguments
    options -- List of command line arguments to parse
    args -- Command line arguments supplied
    """

    #Extract flags from options
    flags = [opt[:-1] for opt in options if opt[-1] == '?']

    try:
        #Add postfix '=' for options that require an argument & add flags without postfix
        long_options = [opt + '=' for opt in options if opt[-1] != '?']
        long_options += flags

        #Call getopt with long arguments only
        tuples, remainder = getopt.getopt(args, '', long_options)
        #If there's a remainder, not all arguments were recognized
        if remainder:
            print '\n'.join(args)
            raise getopt.GetoptError('Unrecognized argument(s) passed: ' + str(remainder), remainder)
        arguments = dict((opt[2:], value) for opt, value in tuples)
    except getopt.GetoptError as err:
        #Print error & usage information to stderr
        print >> sys.stderr, str(err)
        print >> sys.stderr, usage
        sys.exit(1)

    #Remove postfixes '=?' and '?' from options, and '=' postfix from flags
    options = [opt[:-2] if opt[-2:] == '=?' else opt[:-1] if opt[-1] == '?' else opt for opt in options]
    flags = [flag[:-1] if flag[-1] == '=' else flag for flag in flags]

    #Correctly set True/False values for flags, regardless of whether flag was already passed as argument or not
    for flag in flags:
        if flag in arguments:
            #Only overwrite with True if value is empty, as optional arguments (flags) can have values as well
            if not arguments[flag]:
                arguments[flag] = True
        else:
            arguments[flag] = False

    #Ensure all arguments were provided
    for opt in options:
        if opt not in arguments:
            print >> sys.stderr, 'Mandatory argument {0} not provided'.format(opt)
            print >> sys.stderr, usage
            sys.exit(1)

    #Retrieve & return file paths from dictionary in order of options
    return [arguments[option] for option in options]


def get_most_recent_gene_name(genomes, sequence_records):
    """Return gene name annotation for most recently updated genome from sequence records in ortholog."""
    ortholog_products = {}
    for record in sequence_records:
        #Sample header line: >58191|NC_010067.1|YP_001569097.1|COG4948MR|some gene name
        values = record.id.split('|')
        genome = values[0]
        gene_name = values[4]
        ortholog_products[genome] = gene_name

    #If there's only a single gene name, look no further
    if len(set(ortholog_products.values())) == 1:
        return ortholog_products.values()[0]

    #When no genomes are found, for instance when all genomes are external genomes, just return the first annotation
    if not genomes:
        return ortholog_products.values()[0]

    #Determine which genome is the most recent by looking at the modification & release dates of published genomes
    #Starting at the newest genomes, return the first gene name annotation we find
    for genome in sorted(genomes, key=lambda x: x['Modified date'] or x['Released date'], reverse=True):
        if genome['RefSeq project ID'] in ortholog_products:
            return ortholog_products[genome['RefSeq project ID']]
        if genome['Project ID'] in ortholog_products:
            return ortholog_products[genome['Project ID']]

    #Shouldn't really happen, but write this clause anyhow
    logging.warn('Could not retrieve gene name annotation based on date; returning first gene name annotation instead')
    return ortholog_products.values()[0]


def find_cogs_in_sequence_records(sequence_records, include_none=False):
    """Find unique COG annotations assigned to sequences within a single alignment."""
    cogs = set()
    for record in sequence_records:
        #Sample header line: >58191|NC_010067.1|YP_001569097.1|COG4948MR|core
        cog = record.id.split('|')[3]
        if cog in cogs:
            continue
        if cog == 'None':
            cog = None
        if cog != None or include_none:
            cogs.add(cog)
    return cogs
