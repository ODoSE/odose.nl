#!/usr/bin/env python
"""Module to create concatemer per genome of orthologs, create a phylogenetic tree and deduce taxa from that tree."""
from Bio import AlignIO, Phylo, SeqIO
from divergence import create_directory, parse_options, extract_archive_of_files
from divergence.select_taxa import select_genomes_by_ids
from divergence.versions import DNADIST, NEIGHBOR
from subprocess import Popen, PIPE, STDOUT
import logging as log
import os.path
import shutil
import sys
import tempfile

def concatemer_per_genome(run_dir, trimmed_sicos):
    """Create a concatemer DNA file per genome containing all aligned & trimmed SICO genes."""
    concatemer_dir = create_directory('concatemers', inside_dir = run_dir)
    log.info('Creating concatemers from {0} SICOs'.format(len(trimmed_sicos)))

    #Open trimmed concatemer write handles
    concatemer_files = []
    write_handles = {}

    #Loop over trimmed sico files to append each sequence to the right concatemer
    for trimmed_sico in trimmed_sicos:
        for seqr in SeqIO.parse(trimmed_sico, 'fasta'):
            #Sample header line: >58191|NC_010067.1|YP_001569097.1|COG4948MR|core                
            project_id = seqr.id.split('|')[0]

            #Try to retrieve write handle from dictionary of cached write handles per genome
            write_handle = write_handles.get(project_id)

            #If not found, create & store write handle on demand
            if not write_handle:
                #Build up output file path for trimmed SICO genes concatemer per genome
                concatemer_file = os.path.join(concatemer_dir, project_id + '.concat.fna')
                concatemer_files.append(concatemer_file)

                #Open write handle
                write_handle = open(concatemer_file, mode = 'w')
                write_handles[project_id] = write_handle

                #Write initial fasta header
                write_handle.write('> {0}|trimmed concatemer\n'.format(project_id))

            #Write / Append sequence for ortholog to genome concatemer
            write_handle.write('{0}\n'.format(str(seqr.seq)))

    #Close genomes trimmed concatemer write handles 
    for write_handle in write_handles.values():
        write_handle.close()

    log.info('Created {0} concatemers from {1} SICOs'.format(len(concatemer_files), len(trimmed_sicos)))

    return sorted(concatemer_files)

def create_super_concatemer(concatemer_files, destination_path):
    """Concatenate individual genome concatemers into a single super-concatemer for easy import into MEGA viewer."""
    with open(destination_path, mode = 'w') as write_handle:
        for concatemer in concatemer_files:
            seqr = SeqIO.read(concatemer, 'fasta')
            SeqIO.write(seqr, write_handle, 'fasta')

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
    process = Popen(DNADIST, cwd = dnadist_dir, stdin = PIPE, stdout = PIPE, stderr = STDOUT)
    process.communicate(input = 'Y\n')

    #Retrieve outputfile
    outfile = os.path.join(dnadist_dir, 'outfile')
    assert os.path.exists(outfile) and 0 < os.path.getsize(outfile), outfile + ' should exist with some content now'
    return outfile

def _run_neighbor(run_dir, distance_file):
    """Run neighbor to generate a tree of the distances in the distance file, and return the generated tree file."""
    neighbor_dir = create_directory('neighbor', inside_dir = run_dir)

    #Copy outfile from dnadist to infile inside neighbor_dir
    shutil.copy(distance_file, os.path.join(neighbor_dir, 'infile'))

    #Actually run neighbor
    process = Popen(NEIGHBOR, cwd = neighbor_dir, stdin = PIPE, stdout = PIPE, stderr = STDOUT)
    process.communicate(input = 'N\nY\n')

    #Retrieve newick tree file
    treefile = os.path.join(neighbor_dir, 'outtree')
    assert os.path.exists(treefile) and 0 < os.path.getsize(treefile), treefile + ' should exist with some content now'
    return treefile

def _read_taxa_from_tree(tree_file):
    """Read tree_file in Newick format to identify the first two clades that split up this tree and their leafs."""
    #Parse tree using BioPython, which interprets the GenBank Project IDs as confidence scores, but that'll do for now.
    phylo_tree = Phylo.read(tree_file, 'newick')

    #Of the full tree retrieve the clades from the root clade, expecting exactly two distinct clades after UPGMA
    clades = phylo_tree.clade.clades
    assert len(clades) == 2, 'Expected two clades as child of tree\'s first clade, but was {0}'.format(len(clades))

    #Get all the leafs for the above two clades in a similar format to the genome_ids
    clade_one = sorted(str(int(leaf.confidence)) for leaf in clades[0].get_terminals())
    clade_two = sorted(str(int(leaf.confidence)) for leaf in clades[1].get_terminals())
    return clade_one, clade_two

def visualize_tree(super_tree_file, id_to_name_map, ascii_tree):
    """Visualize the phylogenetic tree encoded in the Newick format super_tree_file, and write graphic to ascii_tree."""
    #Draw phylogenetic tree
    tree = Phylo.read(super_tree_file, 'newick')

    #BioPython misinterprets numerical leaf names as confidence scores: Fix this here
    for leaf in tree.get_terminals():
        project_id = str(int(leaf.confidence))
        leaf.name = '{0}\t{1}'.format(project_id, id_to_name_map[project_id])

    #Print ascii tree, as we can't get visualization to work properly using draw_graphviz
    with open(ascii_tree, mode = 'w') as write_handle:
        Phylo.draw_ascii(tree, file = write_handle, column_width = 120)

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: concatenate_orthologs.py
--orthologs-zip=FILE    archive of orthologous genes in FASTA format
--concatemer=FILE       destination file path for super-concatemer of all genomes
--taxon-a=FILE          destination file path for genome IDs for taxon A
--taxon-b=FILE          destination file path for genome IDs for taxon B
--tree=FILE             destination file path for tree visualization
"""
    options = ['orthologs-zip', 'concatemer', 'taxon-a', 'taxon-b', 'tree']
    orthologs_zip, target_concat_file, target_taxon_a, target_taxon_b, target_tree = parse_options(usage, options, args)

    #Run filtering in a temporary folder, to prevent interference from simultaneous runs
    run_dir = tempfile.mkdtemp(prefix = 'concatemer_tree_')

    #Extract files from zip archive
    temp_dir = create_directory('orthologs', inside_dir = run_dir)
    ortholog_files = extract_archive_of_files(orthologs_zip, temp_dir)

    #Concatenate trimmed_files per genome
    concatemer_files = concatemer_per_genome(run_dir, ortholog_files)
    #Create super concatemer
    create_super_concatemer(concatemer_files, target_concat_file)

    #Determine the taxa present in the super concatemer tree by building a phylogenetic tree from genome concatemer and
    #reading genome ids in the two largest clades.
    super_distance_file = _run_dna_dist(run_dir, target_concat_file)
    super_tree_file = _run_neighbor(run_dir, super_distance_file)
    genome_ids_a, genome_ids_b = _read_taxa_from_tree(super_tree_file)

    #Map Project IDs to Organism names
    id_to_name_map = dict((id, genome['Organism Name'])
                          for id, genome in select_genomes_by_ids(genome_ids_a + genome_ids_b).iteritems())

    #Write Project IDs and Organism Names to files
    with open(target_taxon_a, mode = 'w') as write_handle:
        for genome_id in genome_ids_a:
            write_handle.write('{id}\t{name}\n'.format(id = genome_id, name = id_to_name_map[genome_id]))
    with open(target_taxon_b, mode = 'w') as write_handle:
        for genome_id in genome_ids_b:
            write_handle.write('{id}\t{name}\n'.format(id = genome_id, name = id_to_name_map[genome_id]))

    #Visualize tree
    visualize_tree(super_tree_file, id_to_name_map, target_tree)

    #Remove unused files to free disk space 
    shutil.rmtree(run_dir)

    #Exit after a comforting log message
    log.info('Produced: \n%s\n%s\n%s\n%s', target_concat_file, target_taxon_a, target_taxon_b, target_tree)

if __name__ == '__main__':
    main(sys.argv[1:])

