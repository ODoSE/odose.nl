#!/usr/bin/env python
"""package test"""

from divergence.select_taxa import select_taxa
from subprocess import Popen, PIPE
import os.path

#RefSeq project IDs for first 10 ecoli & salmonella
ECOLI_IDS = tuple(str(nr) for nr in [59245, 58531, 59383, 58623, 58783, 58803, 59391, 57915, 58395, 59379])
SALMO_IDS = tuple(str(nr) for nr in [58191, 59431, 58017, 58917, 59247, 59249, 58973, 58831, 59269, 58201])

#Genomes for above IDs retrieved through select_taxa 
ECOLI_GENOMES, SALMONELLA_GENOMES = select_taxa(ECOLI_IDS, SALMO_IDS)
for genome in ECOLI_GENOMES:
    assert genome['RefSeq project ID'] in ECOLI_IDS, '\n'.join(ECOLI_GENOMES)
for genome in SALMONELLA_GENOMES:
    assert genome['RefSeq project ID'] in SALMO_IDS, '\n'.join(SALMONELLA_GENOMES)
ECOLI_AND_SALMONELLA_GENOMES = ECOLI_GENOMES + SALMONELLA_GENOMES
assert len(ECOLI_AND_SALMONELLA_GENOMES) == 20, '\n'.join(ECOLI_AND_SALMONELLA_GENOMES)

def fail_if_diff(testclass, expected_file, actual_file):
    """Run fail_if_diff to compare two files side-by-side, and return exit code & captured standard output."""
    testclass.failUnless(os.path.isfile(expected_file), 'Expected file does not exist' + expected_file)
    testclass.failUnless(0 < os.path.getsize(expected_file), 'Expected file is empty' + expected_file)
    testclass.failUnless(os.path.isfile(actual_file), 'Actual file does not exist' + actual_file)
    testclass.failUnless(0 < os.path.getsize(actual_file), 'Actual file is empty' + actual_file)

    diff_cmd = ['diff', '--side-by-side', '--left-column', expected_file, actual_file]
    process = Popen(diff_cmd, stdout = PIPE)
    stdout = process.communicate()[0]
    testclass.failIf(process.returncode, ' '.join(diff_cmd) + '\n' + stdout)
