#!/usr/bin/env python
"""Module to keep track of paths and versions of other software used within the workflow at various intervals."""

SOFTWARE_DIR = '/projects/divergence/software/'

#Blast
NCBI_BLAST_DIR = SOFTWARE_DIR + 'ncbi-blast-2.2.24+/bin/'
MAKEBLASTDB = NCBI_BLAST_DIR + 'makeblastdb'
BLASTP = NCBI_BLAST_DIR + 'blastp'
BLASTN = NCBI_BLAST_DIR + 'blastn'

#OrthoMCL
ORTHOMCL_DIR = SOFTWARE_DIR + 'orthomclSoftware-v2.0.2/bin/'
ORTHOMCL_INSTALL_SCHEMA = ORTHOMCL_DIR + 'orthomclInstallSchema'
ORTHOMCL_ADJUST_FASTA = ORTHOMCL_DIR + 'orthomclAdjustFasta'
ORTHOMCL_FILTER_FASTA = ORTHOMCL_DIR + 'orthomclFilterFasta'
ORTHOMCL_BLAST_PARSER = ORTHOMCL_DIR + 'orthomclBlastParser'
ORTHOMCL_LOAD_BLAST = ORTHOMCL_DIR + 'orthomclLoadBlast'
ORTHOMCL_PAIRS = ORTHOMCL_DIR + 'orthomclPairs'
ORTHOMCL_DUMP_PAIRS_FILES = ORTHOMCL_DIR + 'orthomclDumpPairsFiles'
ORTHOMCL_MCL_TO_GROUPS = ORTHOMCL_DIR + 'orthomclMclToGroups'
MCL = SOFTWARE_DIR + 'mcl-10-201/src/shmcl/mcl'

#Align & Trim
TRANSLATORX = SOFTWARE_DIR + 'translatorx/translatorx_v1.1.pl'

#Concatemer tree
DNADIST = SOFTWARE_DIR + 'phylip-3.69/exe/dnadist'
NEIGHBOR = SOFTWARE_DIR + 'phylip-3.69/exe/neighbor'

#Calculation
CODEML = SOFTWARE_DIR + 'paml44/bin/codeml'

if __name__ == '__main__':
    pass
