#!/usr/bin/env python
"""Module to calculate pn ps."""
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def calculate_pnps(genomes_a, genomes_b, sico_files):
    """"""
    refseq_ids_a = [genome['RefSeq project ID'] for genome in genomes_a]
    refseq_ids_b = [genome['RefSeq project ID'] for genome in genomes_b]

    print refseq_ids_a
    print refseq_ids_b

    #Separate genomes_a from genomes_b in sico_files
    alignments = (AlignIO.read(sico_file, 'fasta') for sico_file in sico_files)
    for ali in alignments:
        id_seqr_tuples = [(seqr.id.split('|')[0], seqr) for seqr in ali]
        alignment_a = MultipleSeqAlignment(seqr for id, seqr in id_seqr_tuples if id in refseq_ids_a)
        alignment_b = MultipleSeqAlignment(seqr for id, seqr in id_seqr_tuples if id in refseq_ids_b)

        print 'Ali A:\n', alignment_a
        print 'Ali B:\n', alignment_b

if __name__ == '__main__':
    pass
