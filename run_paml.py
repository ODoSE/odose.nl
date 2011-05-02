#!/usr/bin/env python
"""Module to run Phylogenetic Analysis by Maximum Likelihood (PAML)."""

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from divergence import create_directory
from subprocess import check_call, STDOUT
import logging as log
import os.path

def run_paml(genomes_a, genomes_b, sico_files):
    """Run PAML for representatives of clades A and B in each of the SICO files, to calculate dN/dS."""
    #Trash existing paml dir first
    create_directory('paml', delete_first = True)

    representative_a = genomes_a[0]['RefSeq project ID']
    representative_b = genomes_b[0]['RefSeq project ID']

    log.info('Running PAML for {0} aligned and trimmed SICOs'.format(len(sico_files)))

    for sico_file in sico_files:
        _run_yn00(representative_a, representative_b, sico_file)

    return None

YN00 = '/projects/divergence/software/paml44/bin/yn00'

def _run_yn00(repr_id_a, repr_id_b, sico_file):
    """Run yn00 from PAML for selected sequence records from sico_file, returning main nexus output file."""
    #Find sequences from above chosen clade representatives in each SICO file
    seqr_a = None
    seqr_b = None
    for seqrecord in SeqIO.parse(sico_file, 'fasta', alphabet = generic_dna):
        refseq_id = seqrecord.id.split('|')[0]
        if refseq_id == repr_id_a:
            seqr_a = seqrecord
        elif refseq_id == repr_id_b:
            seqr_b = seqrecord
        if seqr_a and seqr_b:
            break
    else:
        assert False, '{0}: Both seqr_a & seqr_b should have been set\n{1}\n{2}'.format(sico_file, seqr_a, seqr_b)

    #Write the representative sequence records out to file in yn00 compatible format
    sico = os.path.splitext(os.path.split(sico_file)[1])[0]
    yn00_dir = create_directory('paml/' + sico)
    nexus_file = os.path.join(yn00_dir, sico + '.nexus')
    _write_nexus_file(seqr_a, seqr_b, nexus_file)

    #Generate yn00 configuration file
    output_file = os.path.join(yn00_dir, 'yn' + sico)
    config_file = os.path.join(yn00_dir, 'yn00.ctl')
    _write_config_file(nexus_file, output_file, config_file)

    #Run yn00
    command = [YN00, os.path.split(config_file)[1]]
    check_call(command, cwd = yn00_dir, stdout = open('/dev/null', mode = 'w'), stderr = STDOUT)

    assert os.path.isfile(output_file) and os.path.getsize(output_file), 'Expected some content in ' + output_file
    return output_file

def _write_nexus_file(seqr_a, seqr_b, nexus_file):
    """Write representative sequences out to a file in the PAML compatible nexus format."""
    nexus_contents = """
#NEXUS
begin data; 
   dimensions ntax={ntax} nchar={nchar}; 
   format datatype=dna missing=? gap=-; 
matrix 
clade_a  {seqa}
clade_b  {seqb}
;
end;
""".format(ntax = 2, nchar = len(seqr_a), seqa = str(seqr_a.seq), seqb = str(seqr_b.seq))
    with open(nexus_file, mode = 'w') as write_handle:
        #SeqIO.write([seqr_a, seqr_b], write_handle, 'nexus')
        write_handle.write(nexus_contents)

def _write_config_file(nexus_file, output_file, config_file):
    """Write a yn00 configuration file using relative paths to the nexus file and output file."""
    config_contents = """
   seqfile = {0} * sequence data file name
   outfile = {1} * main result file
   verbose = 0  * 1: detailed output (list sequences), 0: concise output

     icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below

 weighting = 0  * weighting pathways between codons (0/1)?
commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)? 
*    ndata = 1

* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
""".format(os.path.split(nexus_file)[1], os.path.split(output_file)[1])
    with open(config_file, mode = 'w') as write_handle:
        write_handle.write(config_contents)
