#!/usr/bin/env python
"""Module to run Phylogenetic Analysis by Maximum Likelihood (codeml)."""

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable
from collections import deque
from divergence import create_directory, extract_archive_of_files, create_archive_of_files, parse_options, \
    CODON_TABLE_ID
from divergence.versions import CODEML
from multiprocessing import Pool
from subprocess import check_call, STDOUT
import logging as log
import os.path
import shutil
import sys
import tempfile

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"

def run_codeml_for_sicos(codeml_dir, genome_ids_a, genome_ids_b, sico_files):
    """Run codeml for representatives of clades A and B in each of the SICO files, to calculate dN/dS."""
    log.info('Running codeml for {0} aligned and trimmed SICOs'.format(len(sico_files)))

    pool = Pool()
    future_files = []
    for sico_file in sico_files:
        #Separate alignments for clade A & clade B genomes 
        ali = AlignIO.read(sico_file, 'fasta')
        alignment_a = MultipleSeqAlignment(seqr for seqr in ali if seqr.id.split('|')[0] in genome_ids_a)
        alignment_b = MultipleSeqAlignment(seqr for seqr in ali if seqr.id.split('|')[0] in genome_ids_b)

        #Create sub directory for this run based on sico_file name
        filename = os.path.split(sico_file)[1]
        #Split off everything starting from the first dot
        base_name = filename[:filename.find('.')]
        sub_dir = create_directory(base_name, inside_dir = codeml_dir)

        #Submit for asynchronous calculation
        ft_codeml_file = pool.apply_async(run_codeml, (sub_dir, alignment_a, alignment_b))
        future_files.append(ft_codeml_file)

    return [ft_codeml_file.get() for ft_codeml_file in future_files]

#Using the standard NCBI Bacterial, Archaeal and Plant Plastid Code translation table (11).
BACTERIAL_CODON_TABLE = CodonTable.unambiguous_dna_by_id.get(CODON_TABLE_ID)

def run_codeml(sub_dir, alignment_a, alignment_b):
    """Run codeml from PAML for selected sequence records from sico_file, returning main nexus output file."""
    #Note on whether or not I should be randomizing the below representative selection:
    #"both alternatives have their advantages - just selecting one strain for the divergence calculation means that you 
    # know exactly which strains the divergence comes from - but if this strain is anomalous then you might get some 
    # strange results. i think i would stick with a single strain" - AEW

    #Select first sequences from each clade as representatives
    ab_alignment = MultipleSeqAlignment([alignment_a[0], alignment_b[0]])

    #Codeml chokes when presented with an sequence containing stopcodons: strip those out
    sequence_a = ''
    sequence_b = ''
    for index in range(0, len(ab_alignment[0]), 3):
        codon_a = str(ab_alignment[0][index:index + 3].seq)
        codon_b = str(ab_alignment[1][index:index + 3].seq)
        if codon_a not in BACTERIAL_CODON_TABLE.stop_codons and codon_b not in BACTERIAL_CODON_TABLE.stop_codons:
            sequence_a += codon_a
            sequence_b += codon_b

    #Write the representative sequence records out to file in codeml compatible format
    base_name = os.path.split(sub_dir)[1]
    nexus_file = os.path.join(sub_dir, base_name + '.nexus')
    _write_nexus_file(sequence_a, sequence_b, nexus_file)

    #Generate codeml configuration file
    output_file = os.path.join(sub_dir, base_name + '.codeml')
    config_file = os.path.join(sub_dir, 'codeml.ctl')
    _write_config_file(nexus_file, output_file, config_file)

    #Run codeml
    command = [CODEML, os.path.split(config_file)[1]]
    check_call(command, cwd = sub_dir, stdout = open('/dev/null', mode = 'w'), stderr = STDOUT)

    assert os.path.isfile(output_file) and os.path.getsize(output_file), 'Expected some content in ' + output_file
    return output_file

def _write_nexus_file(sequence_a, sequence_b, nexus_file):
    """Write representative sequences out to a file in the codeml compatible nexus format."""
    nexus_contents = '''
#NEXUS
begin data; 
   dimensions ntax={ntax} nchar={nchar}; 
   format datatype=dna missing=? gap=-; 
matrix 
clade_a  {seqa}
clade_b  {seqb}
;
end;
'''.format(ntax = 2, nchar = len(sequence_a), seqa = sequence_a, seqb = sequence_b)
    with open(nexus_file, mode = 'w') as write_handle:
        write_handle.write(nexus_contents)

def _write_config_file(nexus_file, output_file, config_file):
    """Write a codeml configuration file using relative paths to the nexus file and output file."""
    config_contents = '''
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
'''.format(os.path.split(nexus_file)[1], os.path.split(output_file)[1])
    with open(config_file, mode = 'w') as write_handle:
        write_handle.write(config_contents)

def parse_codeml_output(codeml_file):
    """Parse last line of codeml output file to read initial values, and calculate Dn & Ds as derived values."""
    with open(codeml_file) as read_handle:
        #Extract & parse last line
        last_line = deque(read_handle).pop()
        #Example lines:
        #t=50.0000  S=    97.9  N=   328.1  dN/dS= 0.0113  dN= 0.7872  dS=69.8724
        #t= 1.0569  S=   387.3  N=   950.7  dN/dS= 0.0236  dN= 0.0272  dS= 1.1503

        iterator = iter(item.strip() for item in last_line.replace('=', ' ').split())
        #Use the same above iterator twice in zip to create pairs from sequential items, which we can feed into dict
        value_dict = dict(zip(iterator, iterator))

        for key, value in value_dict.iteritems():
            value_dict[key] = float(value)

        #Below calculations according to AEW to get large D values
        value_dict['Dn'] = float(value_dict['dN']) * float(value_dict['N'])
        value_dict['Ds'] = float(value_dict['dS']) * float(value_dict['S'])
        return value_dict

def _write_dnds_per_ortholog(dnds_file, codeml_files):
    """For each codeml output file write dN, dS & dN/dS to single tab separated file, each on a new line."""
    #Open file to write dN dS values to
    with open(dnds_file, mode = 'w') as write_handle:
        write_handle.write('#Ortholog\tN\tdN\tDn\tS\tdS\tDs\tdN/dS\n')
        #small d and p stand for numbers per site - dN or dn is the number of non-synonymous substitutions per site
        #and Dn is the total number of non-synonymous substitutions

        #Write on each line: SICO file, N, dN, Dn, S, dS, Ds & dN/dS
        for codeml_file in codeml_files:
            sico = os.path.split(codeml_file)[1].split('.')[0]
            value_dict = parse_codeml_output(codeml_file)
            write_handle.write('{0}\t{1[N]}\t{1[dN]}\t{1[Dn]}\t{1[S]}\t{1[dS]}\t{1[Ds]}\t{1[dN/dS]}\n'
                               .format(sico, value_dict))
    return dnds_file

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: run_codeml.py 
--genomes-a=FILE     file with GenBank Project IDs from complete genomes table on each line for taxon A
--genomes-b=FILE     file with GenBank Project IDs from complete genomes table on each line for taxon B
--sico-zip=FILE      archive of aligned & trimmed single copy orthologous (SICO) genes
--codeml-zip=FILE     destination file path for archive of codeml output per SICO gene
--dnds-stats=FILE     destination file path for file with dN, dS & dN/dS values per SICO gene
"""
    options = ['genomes-a', 'genomes-b', 'sico-zip', 'codeml-zip', 'dnds-stats']
    genome_a_ids_file, genome_b_ids_file, sico_zip, codeml_zip, dnds_file = parse_options(usage, options, args)

    #Parse file to extract GenBank Project IDs
    with open(genome_a_ids_file) as read_handle:
        genome_ids_a = [line.split()[0] for line in read_handle]
    with open(genome_b_ids_file) as read_handle:
        genome_ids_b = [line.split()[0] for line in read_handle]

    #Create run_dir to hold files relating to this run
    run_dir = tempfile.mkdtemp(prefix = 'run_codeml_')

    #Extract files from zip archive
    sico_files = extract_archive_of_files(sico_zip, create_directory('sicos', inside_dir = run_dir))

    #Actually run codeml
    codeml_files = run_codeml_for_sicos(run_dir, genome_ids_a, genome_ids_b, sico_files)

    #Write dnds values to single output file
    _write_dnds_per_ortholog(dnds_file, codeml_files)

    #Write the produced files to command line argument filenames
    create_archive_of_files(codeml_zip, codeml_files)

    #Remove unused files to free disk space 
    shutil.rmtree(run_dir)

    #Exit after a comforting log message
    log.info("Produced: \n%s\n%s", codeml_zip, dnds_file)

if __name__ == '__main__':
    main(sys.argv[1:])
