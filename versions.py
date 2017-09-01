#!/usr/bin/env python
"""Module to keep track of paths and versions of other software used within the workflow at various intervals."""

import Bio
import argparse
import getpass
import logging as log
import os
from subprocess import Popen, PIPE, check_call
import sys


__author__ = "Tim te Beek"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


SOFTWARE_DIR = os.path.dirname(os.path.abspath(__file__)) + '/'
if not os.path.isdir(SOFTWARE_DIR):
    logging.error('Software directory is missing: %s', SOFTWARE_DIR)

# Blast
# ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST
MAKEBLASTDB = SOFTWARE_DIR + 'makeblastdb'
BLASTP = SOFTWARE_DIR + 'blastp'
BLASTN = SOFTWARE_DIR + 'blastn'

# Life Science Grid Portal
# https://apps.grid.sara.nl/applications/makeblastdb/
LSGP_MAKEBLASTDB = SOFTWARE_DIR + 'makeblastdb/2.2.26'
LSGP_BLASTN = SOFTWARE_DIR + 'blastn/2.2.26'
LSGP_BLASTP = SOFTWARE_DIR + 'blastp/2.2.27'

# OrthoMCL
# http://www.orthomcl.org/common/downloads/software/v2.0/
ORTHOMCL_INSTALL_SCHEMA = SOFTWARE_DIR + 'orthomclInstallSchema'
ORTHOMCL_ADJUST_FASTA = SOFTWARE_DIR + 'orthomclAdjustFasta'
ORTHOMCL_FILTER_FASTA = SOFTWARE_DIR + 'orthomclFilterFasta'
ORTHOMCL_BLAST_PARSER = SOFTWARE_DIR + 'orthomclBlastParser'
ORTHOMCL_LOAD_BLAST = SOFTWARE_DIR + 'orthomclLoadBlast'
ORTHOMCL_PAIRS = SOFTWARE_DIR + 'orthomclPairs'
ORTHOMCL_DUMP_PAIRS_FILES = SOFTWARE_DIR + 'orthomclDumpPairsFiles'
ORTHOMCL_MCL_TO_GROUPS = SOFTWARE_DIR + 'orthomclMclToGroups'
# http://micans.org/mcl/
MCL = SOFTWARE_DIR + 'mcl'

# Align & Trim
# http://pc16141.mncn.csic.es/cgi-bin/translatorx_vLocal.pl
TRANSLATORX = SOFTWARE_DIR + 'translatorx'

# Concatemer tree
# http://evolution.genetics.washington.edu/phylip.html
PHYLIP = SOFTWARE_DIR + 'phylip'
DNADIST = PHYLIP + ' ' + 'dnadist'
NEIGHBOR = PHYLIP + ' ' + 'neighbor'

# Recombination
# http://www.maths.otago.ac.nz/~dbryant/software.html
PHIPACK = SOFTWARE_DIR + 'Phi'

# Calculation
# http://abacus.gene.ucl.ac.uk/software/paml.html
#PAML_DIR = SOFTWARE_DIR + 'paml4.7/'
#CODEML = PAML_DIR + 'bin/codeml'
CODEML = SOFTWARE_DIR + 'codeml'

def _call_program(*command):
    """Execute command and return the standard output returned by the program. Standard error is caught and ignored."""
    logging.debug(' '.join(command))
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    process.wait()
    return process.communicate()[0].strip()


def _grep_version(path, pattern='version'):
    """Grep for the pattern `version` case insensitively in files specified on path and return the first line."""
    stdout = _call_program('grep', '-ri', pattern, path)
    return stdout.split('\n')[0]


def _parse_args():
    '''
    Parse required arguments.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('target',
                        help='Target output file for version numbers',
                        type=lambda path: logging.FileHandler(path, mode='w'))
    args = parser.parse_args()

    # Directly configure logging through args
    logging.basicConfig(level=logging.INFO, stream=sys.stdout)
    args.target.setFormatter(logging.Formatter())
    logging.root.addHandler(args.target)

    # Return any other args
    return args

def _check_package(pkg_name):
    command = ['which', pkg_name]

    log.info('Executing: %s', ' '.join(command))
    check_call(command, stdout=None)

def main():
    """Method intended to be run when __name-- == '__main__'."""
    # BioPython
    log.info('BioPython\t%s', Bio.__version__)

    # Blast
    print(MAKEBLASTDB)
    _check_package(MAKEBLASTDB)
    _check_package(BLASTP)
    _check_package(BLASTN)

    # Life Science Grid Portal
    #logging.info('LSGP %s', LSGP_MAKEBLASTDB.replace('/', '\t'))
    #logging.info('LSGP %s', LSGP_BLASTP.replace('/', '\t'))
    #logging.info('LSGP %s', LSGP_BLASTN.replace('/', '\t'))

    # OrthoMCL & mcl
    _check_package(ORTHOMCL_INSTALL_SCHEMA)
    _check_package(ORTHOMCL_ADJUST_FASTA)
    _check_package(ORTHOMCL_FILTER_FASTA)
    _check_package(ORTHOMCL_BLAST_PARSER)
    _check_package(ORTHOMCL_LOAD_BLAST)
    _check_package(ORTHOMCL_PAIRS)
    _check_package(ORTHOMCL_DUMP_PAIRS_FILES)
    _check_package(ORTHOMCL_MCL_TO_GROUPS)
    _check_package(MCL)

    # TranslatorX calls muscle internally
    _check_package(TRANSLATORX)
    #logging.info('Muscle\t%s', _call_program('muscle', '-version'))

    # PHYLIP dnadist & neighbor
    _check_package(PHYLIP)
    _check_package(DNADIST)
    _check_package(NEIGHBOR)

    # PHIPACK
    _check_package(PHIPACK)

    # PAML codeml
    _check_package(CODEML)

if __name__ == '__main__':
    # Parse arguments to setup logging; not in main for testing
    #_parse_args()
    # Log software versions
    main()
