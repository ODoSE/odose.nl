#!/usr/bin/env python
"""Module to keep track of paths and versions of other software used within the workflow at various intervals."""

from subprocess import Popen, PIPE
import divergence  # @UnusedImport # pylint: disable=W0611
import Bio
import getpass
import logging
import os

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"

SOFTWARE_DIR = '/work/odosenl/software/'
if getpass.getuser() == 'tbeek':
    SOFTWARE_DIR = '/data/projects/odosenl/software/'
assert os.path.isdir(SOFTWARE_DIR), 'Software directory is missing'

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
    print ' '.join(command)
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    process.wait()
    return process.communicate()[0]


def _grep_version(path, pattern='version'):
    """Grep for the pattern `version` case insensitively in files specified on path and return the first line."""
    stdout = _call_program('grep', '-ri', pattern, path)
    return stdout.split('\n')[0]


def main():
    """Method intended to be run when __name-- == '__main__'."""
    # BioPython
    logging.info('BioPython: ' + Bio.__version__ + '\n')

    # Blast
    logging.info(_call_program(MAKEBLASTDB, '-version'))
    logging.info(_call_program(BLASTP, '-version'))
    logging.info(_call_program(BLASTN, '-version'))

    # Life Science Grid Portal
    logging.info('Life Science Grid Portal BLAST applications: ' + LSGP_BLAST_VERSION)

    # OrthoMCL & mcl
    logging.info('OrthoMCL: ' + _grep_version(ORTHOMCL_DIR + '../doc/OrthoMCLEngine/Main/UserGuide.txt') + '\n')
    logging.info(_call_program(MCL, '--version'))

    # PAML codeml
    logging.info('PAML codeml: ' + _grep_version(PAML_DIR + 'src/paml.h') + '\n')

    # PHYLIP dnadist & neighbor
    logging.info('PHYLIP dnadist: ' + _grep_version(PHYLIP_DIR + 'src/dnadist.c')[3:] + '\n')
    logging.info('PHYLIP neighbor: ' + _grep_version(PHYLIP_DIR + 'src/neighbor.c')[3:] + '\n')

    # TranslatorX calls muscle internally
    logging.info('translatorx: ' + _grep_version(TRANSLATORX, pattern='TranslatorX v')[28:-6] + '\n')
    logging.info(_call_program('muscle', '-version'))

if __name__ == '__main__':
    main()
