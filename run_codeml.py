#!/usr/bin/env python
"""Module to run Phylogenetic Analysis by Maximum Likelihood (codeml)."""

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from divergence import create_directory, extract_archive_of_files, create_archive_of_files, parse_options
from subprocess import check_call, STDOUT
import logging as log
import os.path
import shutil
import sys
import tempfile

def run_codeml(codeml_dir, genome_ids_a, genome_ids_b, sico_files):
    """Run codeml for representatives of clades A and B in each of the SICO files, to calculate dN/dS."""
    #Pick the first genomes as representatives for each clade
    representative_a = genome_ids_a[0]
    representative_b = genome_ids_b[0]

    log.info('Running codeml for {0} aligned and trimmed SICOs'.format(len(sico_files)))
    return [_run_codeml(codeml_dir, representative_a, representative_b, sico_file) for sico_file in sico_files]

CODEML = '/projects/divergence/software/paml44/bin/codeml'

def _run_codeml(codeml_dir, repr_id_a, repr_id_b, sico_file):
    """Run codeml from PAML for selected sequence records from sico_file, returning main nexus output file."""
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

    #Write the representative sequence records out to file in codeml compatible format
    sico = os.path.splitext(os.path.split(sico_file)[1])[0]
    codeml_dir = create_directory(sico, inside_dir = codeml_dir)
    nexus_file = os.path.join(codeml_dir, sico + '.nexus')
    _write_nexus_file(seqr_a, seqr_b, nexus_file)

    #Generate codeml configuration file
    output_file = os.path.join(codeml_dir, 'codeml_' + sico)
    config_file = os.path.join(codeml_dir, 'codeml.ctl')
    _write_config_file(nexus_file, output_file, config_file)

    #Run codeml
    command = [CODEML, os.path.split(config_file)[1]]
    check_call(command, cwd = codeml_dir, stdout = open('/dev/null', mode = 'w'), stderr = STDOUT)

    assert os.path.isfile(output_file) and os.path.getsize(output_file), 'Expected some content in ' + output_file
    return output_file

def _write_nexus_file(seqr_a, seqr_b, nexus_file):
    """Write representative sequences out to a file in the codeml compatible nexus format."""
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
    """Write a codeml configuration file using relative paths to the nexus file and output file."""
    config_contents = """
      seqfile = {0} * sequence data filename
      outfile = {1}           * main result file name
     treefile = test.tree      * tree structure file name

        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0  * 1: detailed output, 0: concise output
      runmode = -2  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table


       clock = 0   * 0:no clock, 1:global clock; 2:local clock; 3:TipDate


       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
   aaRatefile = wag.dat * only used for aa seqs with model=empirical(_F)
                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

        model = 0
                   * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                   * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
        Mgene = 1  * 0:rates, 1:separate; 

    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2  * initial or fixed kappa
    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate 
        omega = .4 * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 8  * # of categories in dG of NSsites models

      fix_rho = 1
          rho = 0.

        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = .5e-6
*   cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
* fix_blength = 0
       method = 0   * 0: simultaneous; 1: one branch at a time
""".format(os.path.split(nexus_file)[1], os.path.split(output_file)[1])
    with open(config_file, mode = 'w') as write_handle:
        write_handle.write(config_contents)

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: run_codeml.py 
--genomes-a=FILE     file with RefSeq id from complete genomes table on each line for clade A
--genomes-b=FILE     file with RefSeq id from complete genomes table on each line for clade B
--sico-zip=FILE      archive of aligned & trimmed single copy orthologous (SICO) genes
--codeml-zip=FILE    destination file path for archive of codeml output per SICO gene
"""
    options = ['genomes-a', 'genomes-b', 'sico-zip', 'codeml-zip']
    genome_a_ids_file, genome_b_ids_file, sico_zip, codeml_zip = parse_options(usage, options, args)

    #Parse file containing RefSeq project IDs to extract RefSeq project IDs
    with open(genome_a_ids_file) as read_handle:
        genome_ids_a = [line.split()[0] for line in read_handle]
    with open(genome_b_ids_file) as read_handle:
        genome_ids_b = [line.split()[0] for line in read_handle]

    #Create run_dir to hold files relating to this run
    run_dir = tempfile.mkdtemp(prefix = 'run_codeml_')

    #Extract files from zip archive
    sico_files = extract_archive_of_files(sico_zip, create_directory('sicos', inside_dir = run_dir))

    #Actually run codeml
    codeml_files = run_codeml(run_dir, genome_ids_a, genome_ids_b, sico_files)

    #Write the produced files to command line argument filenames
    create_archive_of_files(codeml_zip, codeml_files)

    #Remove unused files to free disk space 
    shutil.rmtree(run_dir)

    #Exit after a comforting log message
    log.info("Produced: \n%s", codeml_zip)
    return codeml_zip

if __name__ == '__main__':
    main(sys.argv[1:])
