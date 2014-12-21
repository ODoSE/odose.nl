#!/usr/bin/env python
"""Module to keep track of paths and versions of other software used within the workflow at various intervals."""

import Bio
import argparse
import getpass
import logging
import os
from subprocess import Popen, PIPE
import sys


__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


SOFTWARE_DIR = '/work/odosenl/software/'
if getpass.getuser() == 'tim':
    SOFTWARE_DIR = '/home/tim/Documents/dev/odosenl/software/'
if not os.path.isdir(SOFTWARE_DIR):
    logging.error('Software directory is missing: %s', SOFTWARE_DIR)

# Blast
NCBI_BLAST_DIR = SOFTWARE_DIR + 'ncbi-blast-2.2.28+/bin/'
MAKEBLASTDB = NCBI_BLAST_DIR + 'makeblastdb'
BLASTP = NCBI_BLAST_DIR + 'blastp'
BLASTN = NCBI_BLAST_DIR + 'blastn'

# Life Science Grid Portal
LSGP_BLAST_VERSION = '2.2.26'
LSGP_MAKEBLASTDB = 'makeblastdb/' + LSGP_BLAST_VERSION
LSGP_BLASTN = 'blastn/' + LSGP_BLAST_VERSION
LSGP_BLASTP = 'blastp/' + LSGP_BLAST_VERSION

# OrthoMCL
MCL = SOFTWARE_DIR + 'mcl-12-135/src/shmcl/mcl'
ORTHOMCL_DIR = SOFTWARE_DIR + 'orthomclSoftware-v2.0.7/bin/'
ORTHOMCL_INSTALL_SCHEMA = ORTHOMCL_DIR + 'orthomclInstallSchema'
ORTHOMCL_ADJUST_FASTA = ORTHOMCL_DIR + 'orthomclAdjustFasta'
ORTHOMCL_FILTER_FASTA = ORTHOMCL_DIR + 'orthomclFilterFasta'
ORTHOMCL_BLAST_PARSER = ORTHOMCL_DIR + 'orthomclBlastParser'
ORTHOMCL_LOAD_BLAST = ORTHOMCL_DIR + 'orthomclLoadBlast'
ORTHOMCL_PAIRS = ORTHOMCL_DIR + 'orthomclPairs'
ORTHOMCL_DUMP_PAIRS_FILES = ORTHOMCL_DIR + 'orthomclDumpPairsFiles'
ORTHOMCL_MCL_TO_GROUPS = ORTHOMCL_DIR + 'orthomclMclToGroups'

# Align & Trim
TRANSLATORX = SOFTWARE_DIR + 'translatorx/translatorx_v1.1.pl'

# Concatemer tree
PHYLIP_DIR = SOFTWARE_DIR + 'phylip-3.69/'
DNADIST = PHYLIP_DIR + 'exe/dnadist'
NEIGHBOR = PHYLIP_DIR + 'exe/neighbor'

# Recombination
PHIPACK = SOFTWARE_DIR + 'PhiPack/Phi'

# Calculation
PAML_DIR = SOFTWARE_DIR + 'paml4.7/'
CODEML = PAML_DIR + 'bin/codeml'


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


def main():
    """Method intended to be run when __name-- == '__main__'."""
    # BioPython
    logging.info('BioPython\t' + Bio.__version__)

    # Blast
    logging.info('makeblastdb\t' + _call_program(MAKEBLASTDB, '-version').split('\n')[1])
    logging.info('makeblastdb\t' + _call_program(BLASTP, '-version').split('\n')[1])
    logging.info('makeblastdb\t' + _call_program(BLASTN, '-version').split('\n')[1])

    # Life Science Grid Portal
    logging.info('LSGP BLAST\t' + LSGP_BLAST_VERSION)

    # OrthoMCL & mcl
    logging.info('OrthoMCL\t' + _grep_version(ORTHOMCL_DIR + '../doc/OrthoMCLEngine/Main/UserGuide.txt'))
    logging.info('MCL\t' + _call_program(MCL, '--version').split('\n')[0])

    # PAML codeml
    logging.info('PAML codeml\t' + _grep_version(PAML_DIR + 'src/paml.h'))

    # PHYLIP dnadist & neighbor
    logging.info('PHYLIP dnadist\t' + _grep_version(PHYLIP_DIR + 'src/dnadist.c')[3:])
    logging.info('PHYLIP neighbor\t' + _grep_version(PHYLIP_DIR + 'src/neighbor.c')[3:])

    # TranslatorX calls muscle internally
    logging.info('translatorx\t' + _grep_version(TRANSLATORX, pattern='TranslatorX v')[28:-6])
    logging.info('Muscle\t' + _call_program('muscle', '-version'))

if __name__ == '__main__':
    # Parse arguments to setup logging; not in main for testing
    _parse_args()
    # Log software versions
    main()
