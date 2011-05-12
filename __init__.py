#!/usr/bin/env python
"""Package divergence"""

from pkg_resources import resource_filename #IGNORE:E0611 @UnresolvedImport
from zipfile import ZipFile, ZIP_DEFLATED
import Bio
import getopt
import httplib2
import logging as log
import os
import shutil
import sys

#Setup basic logging configuration with log level INFO  
log.basicConfig(level = log.INFO,
                stream = sys.stdout,
                format = '%(levelname)s\t%(asctime)s %(module)s.%(funcName)s:%(lineno)d\t%(message)s',
                datefmt = '%H:%M:%S')

#Require at least version 1.53 op BioPython
assert 1.54 <= float(Bio.__version__), 'BioPython version 1.54 or higher is required'

#Base output dir
BASE_OUTPUT_PATH = '../divergence-cache/'

def create_directory(dirname, inside_dir = BASE_OUTPUT_PATH):
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
    
    WARNING: The binary mode implies new \n characters will NOT be added in between files!"""
    with open(target_path, mode = 'wb') as write_handle:
        for source_file in source_files:
            shutil.copyfileobj(open(source_file, mode = 'rb'), write_handle)
    assert os.path.isfile(target_path) and 0 < os.path.getsize(target_path), target_path + ' should exist with content'

def create_archive_of_files(target_archive_file, file_iterable):
    """Write all files contained in file_iterable to target_archive_file, with archived file names truncated to last part."""
    write_handle = ZipFile(target_archive_file, mode = 'w', compression = ZIP_DEFLATED)
    for some_file in file_iterable:
        write_handle.write(some_file, os.path.split(some_file)[1])
    write_handle.close()

def extract_archive_of_files(archive_file, target_dir):
    """Extract all files from archive to target directory, and return list of files extracted."""
    extracted_files = []
    read_handle = ZipFile(archive_file, mode = 'r')
    for zipinfo in read_handle.infolist():
        target_path = os.path.join(target_dir, zipinfo.filename)
        extracted_path = read_handle.extract(zipinfo, path = target_path)
        extracted_files.append(extracted_path)
    read_handle.close()
    return extracted_files

def parse_options(usage, options, args):
    """Parse command line arguments in args, using mandatory options  
    
    Arguments:
    usage -- Usage string detailing command line arguments
    options -- List of command line arguments to parse
    args -- Command line arguments supplied
    """
    try:
        #postfix '=' to indicate options require an argument
        long_options = [opt + '=' for opt in options]
        tuples, remainder = getopt.getopt(args, '', long_options)
        if remainder:
            raise getopt.GetoptError('Unrecognized argument(s) passed', remainder)
        arguments = dict((opt[2:], value) for opt, value in tuples)
    except getopt.GetoptError as err:
        #Print error & usage information to stderr 
        print >> sys.stderr, str(err)
        print >> sys.stderr, usage
        sys.exit(1)

    #Ensure all arguments were provided
    for opt in options:
        if opt not in arguments:
            print >> sys.stderr, 'Mandatory argument {0} not provided'.format(opt)
            print >> sys.stderr, usage
            sys.exit(1)

    #Retrieve & return file paths from dictionary
    return [arguments[option] for option in options]

