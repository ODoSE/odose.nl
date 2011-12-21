#!/usr/bin/env python
"""Module to filter orthologs either with multiple COG annotations or when recombination is found."""

from __future__ import division
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from divergence import create_directory, extract_archive_of_files, create_archive_of_files, parse_options
from divergence.concatemer_tree import _run_dna_dist, _run_neighbor, _read_taxa_from_tree
from divergence.crosstable_gene_ids import create_crosstable
import logging as log
import os.path
import shutil
import sys
import tempfile

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


def _filter_multiple_cog_orthologs(run_dir, ortholog_files):
    """Filter orthologs where multiple different COG annotations are found, and in addition transfer COGs."""

    log.info('Filtering orthologs with multiple COG annotations')

    #Retrieve SICO to cog dictionaries of cog conflicts & transferable cog annotations and list of SICOs missing cog
    cog_conflicts, cog_transferable, cog_missing = _group_cog_issues(ortholog_files)
    _log_cog_statistics(cog_conflicts, cog_transferable, cog_missing)

    #Filter out orthologs containing more than one COG annotation
    ortholog_files = [sico for sico in ortholog_files if sico not in cog_conflicts.keys()]

    #File detailing transfered COG annotations for recipient protein IDs & COGs
    transfered_cogs = os.path.join(run_dir, 'transfered_cogs.tsv')
    with open(transfered_cogs, mode='w') as write_handle:
        write_handle.write('ProjectID\tAccessioncode\tProteinID\tCOG\tsource\n')

    #Transfer COGs by overwriting sico_files with correct COG set
    for sico_file, cog in cog_transferable.iteritems():
        _assign_cog_to_sequences(sico_file, cog, transfered_cogs)

    return ortholog_files, transfered_cogs


def find_cogs_in_sequence_records(sequence_records, include_none=False):
    """Find unique COG annotations assigned to sequences within a single alignment."""
    cogs = set()
    for record in sequence_records:
        #Sample header line: >58191|NC_010067.1|YP_001569097.1|COG4948MR|core
        cog = record.id.split('|')[3]
        if cog in cogs:
            continue
        if cog == 'None':
            cog = None
        if cog != None or include_none:
            cogs.add(cog)
    return cogs


def _group_cog_issues(sico_files):
    """Find issues with COG assignments within SICO files by looking at COG conflicts, transferable and missing COGs."""
    cog_conflicts = {}
    cog_transferable = {}
    cog_missing = []
    for sico_file in sico_files:
        cogs = find_cogs_in_sequence_records(SeqIO.parse(sico_file, 'fasta'), include_none=True)
        if 0 == len(cogs):
            cog_missing.append(sico_file)
            continue
        if 1 < len(cogs):
            if None in cogs:
                cogs.remove(None)
                if len(cogs) == 1:
                    cog_transferable[sico_file] = cogs.pop()
                    continue
            cog_conflicts[sico_file] = cogs
    return cog_conflicts, cog_transferable, cog_missing


def _assign_cog_to_sequences(fasta_file, cog, transfered_cogs):
    """Assign cog annotatione to all sequences in fasta_file."""
    #Read all sequence records from file
    seqrecords = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta')).values()
    #Write all sequence records back to the same(!) file
    with open(fasta_file, mode='w') as write_handle:
        with open(transfered_cogs, mode='a') as append_handle:
            for seqr in seqrecords:
                #Sample header line: >58191|NC_010067.1|YP_001569097.1|COG4948MR|core
                #Or for missing COG: >58191|NC_010067.1|YP_001569097.1|None|core
                split = seqr.id.split('|')
                assert split[3] in (cog, 'None'), 'COG should be either {0} or None, but was {1}'.format(cog, split[3])
                if split[3] == 'None':
                    #Assign cog and alter seqr variable to include assigned cog
                    split[3] = cog
                    seqr = SeqRecord(seqr.seq, id='|'.join(split), description='')
                    #Append this COG transfer to append_handle
                    append_handle.write('\t'.join(split) + '\n')
                SeqIO.write(seqr, write_handle, 'fasta')


def _log_cog_statistics(cog_conflicts, cog_transferable, cog_missing):
    """Append COG statistics to stats_file"""
    if cog_conflicts:
        log.info('{0}\tOrthologs contained multiple COGs and were therefore removed:'.format(len(cog_conflicts)))
        for sico_file in sorted(cog_conflicts.keys()):
            cogs = cog_conflicts[sico_file]
            log.info('{0}:\t{1}'.format(os.path.split(sico_file)[1], '\t'.join(cogs)))
    if cog_transferable:
        log.info('{0}\tCOGs transfered within orthologs to unannotated sequences'.format(len(cog_transferable)))
    if cog_missing:
        log.info('{0}\tOrthologs did not contain any COG annotations'.format(len(cog_missing)))


def _filter_recombined_orthologs(run_dir, aligned_files, genome_ids_a, genome_ids_b):
    """Filter aligned fasta files where there is evidence of recombination when inspecting phylogenetic trees.
    Return two collections of aligned files, the first without recombination, the second with recombination."""

    log.info('Filtering orthologs where phylogenetic trees show evidence of recombination')

    #Collections to hold both non recombination files & files showing recombination
    non_recomb = []
    recombined = []

    #Assign ortholog files to the correct collection based on whether they show recombination
    for ortholog_file in aligned_files:
        #Determine input file base name to create an ortholog run specific directory
        base_name = os.path.split(os.path.splitext(ortholog_file)[0])[1]
        ortholog_dir = create_directory(base_name, inside_dir=run_dir)

        #Create distance file
        distance_file = _run_dna_dist(ortholog_dir, ortholog_file)

        #Create tree file
        tree_file = _run_neighbor(ortholog_dir, distance_file)

        #Parse tree file to ensure all genome_ids_a & genome_ids_b group together in the tree
        if _tree_shows_recombination(genome_ids_a, genome_ids_b, tree_file):
            recombined.append(ortholog_file)
        else:
            non_recomb.append(ortholog_file)

    log.info('%i Orthologs out of %i were filtered out due to recombination, leaving %i non recombined orthologs',
             len(recombined), len(aligned_files), len(non_recomb))

    return non_recomb, recombined


def _tree_shows_recombination(genome_ids_a, genome_ids_b, tree_file):
    """Look for evidence of recombination by seeing if all genomes of the separate taxa group together in the tree."""

    #Sample tree: (((((59245:0.00000,58803:0.00000):0.00000,59391:0.00000):0.01222,
    #((58531:0.00000,57915:0.00000):0.00182,(58623:0.00091,59379:0.00091):0.00091):0.01040):0.00229,
    #((59383:0.00457,58395:0.00457):0.00601,58783:0.01059):0.00393):0.06715,(58191:0.02967,(((59431:0.00367,
    #(58973:0.00000,58831:0.00000):0.00367):0.00117,(((58917:0.00000,59247:0.00000):0.00000,59249:0.00000):0.00367,
    #(59269:0.00000,58201:0.00000):0.00367):0.00117):0.00172,58017:0.00655):0.02311):0.05199)

    clade_one, clade_two = _read_taxa_from_tree(tree_file)

    #Use first genome of clade A to determine which collections should match with one another
    first_a_id = genome_ids_a[0]
    if first_a_id in clade_one:
        #We'll declare to have found recombination when the taxa identified through the tree do not match the user taxa
        return set(genome_ids_a) != set(clade_one) or set(genome_ids_b) != set(clade_two)

    assert first_a_id in clade_two, '{0}\n{1}\n{2}\n{3}'.format(tree_file, clade_one, clade_two, first_a_id)
    return set(genome_ids_a) != set(clade_two) or set(genome_ids_b) != set(clade_one)


def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: filter_orthologs.py
--orthologs-zip=FILE            archive of orthologous genes in FASTA format
--filter-multiple-cogs          filter orthologs with multiple COG annotations among genes [OPTIONAL]
--filter-recombination=FILE     filter orthologs that show recombination when comparing phylogenetic trees [OPTIONAL]
                                destination file path for archive of recombination orthologs
--recombined-crosstable=FILE    destination file path for recombined crosstable of GeneIDs, COGs and Products [OPTIONAL]
--taxon-a=FILE                  file with genome IDs for taxon A to use in recombination filtering
--taxon-b=FILE                  file with genome IDs for taxon B to use in recombination filtering
--retained-zip=FILE             destination file path for archive of retained orthologs after filtering
"""
    options = ('orthologs-zip', 'filter-multiple-cogs=?', 'filter-recombination=?', 'recombined-crosstable=?',
               'taxon-a=?', 'taxon-b=?', 'retained-zip')
    orthologs_zip, filter_cogs, filter_recombination, recombined_crosstable, \
    taxona, taxonb, retained_zip = parse_options(usage, options, args)

    #Run filtering in a temporary folder, to prevent interference from simultaneous runs
    run_dir = tempfile.mkdtemp(prefix='filter_orthologs_')

    #Extract files from zip archive
    temp_dir = create_directory('orthologs', inside_dir=run_dir)
    ortholog_files = extract_archive_of_files(orthologs_zip, temp_dir)

    #Filter orthologs with multiple COG annotations among genes if flag was set
    if filter_cogs:
        ortholog_files, transfered_cogs = _filter_multiple_cog_orthologs(run_dir, ortholog_files)

    #Possible extension: filter ortholog when any strain has been flagged as 'mobile element', 'phage' or 'IS element'

    #Filter orthologs that show recombination when comparing phylogenetic trees if flag was set
    if filter_recombination:
        #Parse file to extract GenBank Project IDs
        with open(taxona) as read_handle:
            genome_ids_a = [line.split()[0] for line in read_handle]
        with open(taxonb) as read_handle:
            genome_ids_b = [line.split()[0] for line in read_handle]
        ortholog_files, recombined_files = _filter_recombined_orthologs(run_dir, ortholog_files,
                                                                       genome_ids_a, genome_ids_b)
        #Create crosstable
        create_crosstable(recombined_files, recombined_crosstable)

    #Create archives of files on command line specified output paths
    if filter_cogs:
        shutil.move(transfered_cogs, filter_cogs)
    if filter_recombination:
        create_archive_of_files(filter_recombination, recombined_files)
    create_archive_of_files(retained_zip, ortholog_files)

    #Remove unused files to free disk space
    shutil.rmtree(run_dir)

    #Exit after a comforting log message
    log.info('Produced:')
    if filter_cogs:
        log.info(filter_cogs)
    if filter_recombination:
        log.info(filter_recombination)
    log.info(retained_zip)

if __name__ == '__main__':
    main(sys.argv[1:])
