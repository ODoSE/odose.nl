#!/usr/bin/env python
###
# Part of the Adaptive Divergence through Direction of Selection workflow.
# Copyright (C) 2011  Tim te Beek <tim.te.beek@nbic.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###
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
