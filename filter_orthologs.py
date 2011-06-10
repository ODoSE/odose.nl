#!/usr/bin/env python
"""Module to filter orthologs either with multiple COG annotations or when recombination is found."""

from __future__ import division
from Bio import AlignIO, Phylo, SeqIO
from Bio.SeqRecord import SeqRecord
from divergence import create_directory, extract_archive_of_files, create_archive_of_files, \
    parse_options
from divergence.concatenate_orthologs import concatemer_per_genome, create_super_concatemer
from subprocess import Popen, PIPE, STDOUT
import logging as log
import os.path
import shutil
import sys
import tempfile

def _filter_multiple_cog_orthologs(ortholog_files):
    """Filter orthologs where multiple different COG annotations are found, and in addition transfer COGs."""

    log.info('Filtering orthologs with multiple COG annotations')

    #Retrieve SICO to cog dictionaries of cog conflicts & transferable cog annotations and list of SICOs missing cog
    cog_conflicts, cog_transferable, cog_missing = _group_cog_issues(ortholog_files)
    _log_cog_statistics(cog_conflicts, cog_transferable, cog_missing)

    #Filter out orthologs containing more than one COG annotation
    ortholog_files = [sico for sico in ortholog_files if sico not in cog_conflicts.keys()]

    #Transfer COGs by overwriting sico_files with correct COG set
    for sico_file, cog in cog_transferable.iteritems():
        #TODO Create file listing transfered COG annotations with donor and recipient protein IDs     
        _assign_cog_to_sequences(sico_file, cog)

    return ortholog_files

def find_cogs_in_sequence_records(sequence_records, include_none = False):
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
        cogs = find_cogs_in_sequence_records(SeqIO.parse(sico_file, 'fasta'), include_none = True)
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

def _assign_cog_to_sequences(fasta_file, cog):
    """Assign cog annotatione to all sequences in fasta_file."""
    seqrecords = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta')).values()
    with open(fasta_file, mode = 'w') as write_handle:
        for seqr in seqrecords: #Sample header line: >58191|NC_010067.1|YP_001569097.1|COG4948MR|core
            #Or for missing COG: >58191|NC_010067.1|YP_001569097.1|None|core
            split = seqr.id.split('|')
            assert split[3] in (cog, 'None'), 'COG should be either {0} or None, but was {1}'.format(cog, split[3])
            split[3] = cog
            seqr = SeqRecord(seqr.seq, id = '|'.join(split), description = '')
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

def filter_recombined_orthologs(run_dir, aligned_files):
    """Filter aligned fasta files where there is evidence of recombination when inspecting phylogenetic trees. 
    Return two collections of aligned files, the first without recombination, the second with recombination."""

    log.info('Filtering orthologs where phylogenetic trees show evidence of recombination')

    #Collections to hold both non recombination files & files showing recombination 
    non_recomb = []
    recombined = []

    #Determine the taxa present in the super concatemer tree
    concatemer_files = concatemer_per_genome(run_dir, aligned_files)
    super_concatemer = os.path.join(run_dir, 'super_concatemer.fna')
    create_super_concatemer(concatemer_files, super_concatemer)
    super_distance_file = _run_dna_dist(run_dir, super_concatemer)
    super_tree_file = _run_neighbor(run_dir, super_distance_file)
    genome_ids_a, genome_ids_b = _read_taxa_from_tree(super_tree_file)

    #TODO Ensure the tree that's created here matches the tree that is created from the filtered dataset later
    #Otherwise fail after print >> stderr, 'Unfiltered & filtered tree clustering does not match'

    #Assign ortholog files to the correct collection based on whether they show recombination
    log.info('Recombination found in the following orthologs:')
    for ortholog_file in aligned_files:
        #Determine input file base name to create an ortholog run specific directory
        base_name = os.path.split(os.path.splitext(ortholog_file)[0])[1]
        ortholog_dir = create_directory(base_name, inside_dir = run_dir)

        #Create distance file
        distance_file = _run_dna_dist(ortholog_dir, ortholog_file)

        #Create tree file
        tree_file = _run_neighbor(ortholog_dir, distance_file)

        #Parse tree file to ensure all genome_ids_a & genome_ids_b group together in the tree
        if _find_recombination(genome_ids_a, genome_ids_b, tree_file):
            recombined.append(ortholog_file)
            log.info('%s\t%s', base_name, '\t'.join(find_cogs_in_sequence_records(SeqIO.parse(ortholog_file, 'fasta'))))
        else:
            non_recomb.append(ortholog_file)

    log.info('%i\tOrthologs out of %i were filtered out due to recombination', len(recombined), len(aligned_files))

    return non_recomb, recombined

DNADIST = '/projects/divergence/software/phylip-3.69/exe/dnadist'

def _run_dna_dist(run_dir, aligned_file):
    """Run dnadist to calculate distances between individual strains in a distance matrix, as input for neighbor."""
    #Run calculations inside a directory
    dnadist_dir = create_directory('dnadist/', inside_dir = run_dir)

    #Read alignment file
    alignment = AlignIO.read(aligned_file, 'fasta')

    #Convert alignment in to proper input file for dnadist according to specification
    nr_of_species = len(alignment)
    nr_of_sites = len(alignment[0])
    infile = os.path.join(dnadist_dir, 'infile')
    with open(infile, mode = 'w') as write_handle:
        write_handle.write('   {0}   {1}\n'.format(nr_of_species, nr_of_sites))

        for seq_record in alignment:
            name = seq_record.id.split('|')[0]
            write_handle.write('{0:10}{1}\n'.format(name, seq_record.seq))

    #Actually run the dnadist program in the correct directory, and send input to it for the first prompt
    log.debug('Executing: %s in %s', DNADIST, dnadist_dir)
    process = Popen(DNADIST, cwd = dnadist_dir, stdin = PIPE, stdout = PIPE, stderr = STDOUT)
    process.communicate(input = 'Y\n')

    #Retrieve outputfile
    outfile = os.path.join(dnadist_dir, 'outfile')
    assert os.path.exists(outfile) and 0 < os.path.getsize(outfile), outfile + ' should exist with some content now'
    return outfile

NEIGHBOR = '/projects/divergence/software/phylip-3.69/exe/neighbor'

def _run_neighbor(run_dir, distance_file):
    """Run neighbor to generate a tree of the distances in the distance file, and return the generated tree file."""
    neighbor_dir = create_directory('neighbor', inside_dir = run_dir)

    #Copy outfile from dnadist to infile inside neighbor_dir
    shutil.copy(distance_file, os.path.join(neighbor_dir, 'infile'))

    #Actually run neighbor
    log.debug('Executing: %s in %s', NEIGHBOR, neighbor_dir)
    process = Popen(NEIGHBOR, cwd = neighbor_dir, stdin = PIPE, stdout = PIPE, stderr = STDOUT)
    process.communicate(input = 'N\nY\n')

    #Retrieve newick tree file
    treefile = os.path.join(neighbor_dir, 'outtree')
    assert os.path.exists(treefile) and 0 < os.path.getsize(treefile), treefile + ' should exist with some content now'
    return treefile

def _read_taxa_from_tree(tree_file):
    """Read tree_file in Newick format to identify the first two clades that split up this tree and their leafs."""
    #Parse tree using BioPython, which wrongly interprets the RefSeq IDs as confidence scores, but that'll do for now.
    phylo_tree = Phylo.read(tree_file, 'newick')

    #Of the full tree retrieve the clades from the root clade, expecting exactly two distinct clades after UPGMA
    clades = phylo_tree.clade.clades
    assert len(clades) == 2, 'Expected two clades as child of tree\'s first clade, but was {0}'.format(len(clades))

    #Get all the leafs for the above two clades in a similar format to the genome_ids
    clade_one = sorted(str(int(leaf.confidence)) for leaf in clades[0].get_terminals())
    clade_two = sorted(str(int(leaf.confidence)) for leaf in clades[1].get_terminals())
    return clade_one, clade_two

def _find_recombination(genome_ids_a, genome_ids_b, tree_file):
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
        recombination_found = set(genome_ids_a) != set(clade_one) or set(genome_ids_b) != set(clade_two)
    else:
        assert first_a_id in clade_two, '{0}\n{1}\n{2}\n{3}'.format(tree_file, clade_one, clade_two, first_a_id)
        recombination_found = set(genome_ids_a) != set(clade_two) or set(genome_ids_b) != set(clade_one)

    return recombination_found

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: filter_orthologs.py
--orthologs-zip=FILE           archive of orthologous genes in FASTA format
--filter-multiple-cogs         filter orthologs with multiple COG annotations among genes [OPTIONAL]
--filter-recombination=FILE    filter orthologs that show recombination when comparing phylogenetic trees [OPTIONAL]
                               destination file path for archive of recombination orthologs
--retained-zip=FILE            destination file path for archive of retained orthologs after filtering
"""
    options = ['orthologs-zip', 'filter-multiple-cogs?', 'filter-recombination=?', 'retained-zip']
    orthologs_zip, filter_cogs_enabled, filter_recombination, retained_zip = parse_options(usage, options, args)

    #Run filtering in a temporary folder, to prevent interference from simultaneous runs
    run_dir = tempfile.mkdtemp(prefix = 'filter_orthologs_')

    #Extract files from zip archive
    temp_dir = create_directory('orthologs', inside_dir = run_dir)
    ortholog_files = extract_archive_of_files(orthologs_zip, temp_dir)

    #Filter orthologs with multiple COG annotations among genes if flag was set
    if filter_cogs_enabled:
        ortholog_files = _filter_multiple_cog_orthologs(ortholog_files)

    #TODO Add option to filter out SICOs when any ortholog has been flagged as 'mobile element', 'phage' or 'IS element'

    #Filter orthologs that show recombination when comparing phylogenetic trees if flag was set
    if filter_recombination:
        ortholog_files, recombined_files = filter_recombined_orthologs(run_dir, ortholog_files)

    #Create archives of files on command line specified output paths & move trim_stats_file
    create_archive_of_files(retained_zip, ortholog_files)
    if filter_recombination:
        #FIXME Galaxy fails the run when no recombination files are found, probably because this file is empty
        create_archive_of_files(filter_recombination, recombined_files)

    #Remove unused files to free disk space 
    shutil.rmtree(run_dir)

    #Exit after a comforting log message
    log.info('Produced: ')
    log.info(retained_zip)
    if filter_recombination:
        log.info(filter_recombination)

    if filter_recombination:
        return retained_zip, filter_recombination
    return retained_zip

if __name__ == '__main__':
    main(sys.argv[1:])
