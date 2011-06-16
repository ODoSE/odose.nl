#!/usr/bin/env python
"""Module to calculate pn ps."""

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable
from divergence import parse_options, create_directory, extract_archive_of_files
from divergence.filter_orthologs import find_cogs_in_sequence_records
from divergence.run_codeml import run_codeml, parse_codeml_output
from itertools import product
import logging as log
import os.path
import re
import shutil
import sys
import tempfile

#TODO Rename this module to something more fitting

def _append_sums_and_dos_average(calculations_file, sfs_max_nton, comp_values_list):
    """Append sums over columns 3 through -1, and the mean of the final direction of selection column."""
    sum_comp_values = {}
    dos_list_a = []
    for comp_values in comp_values_list:
    #Sum the following columns
        for column in _get_column_headers_in_sequence(sfs_max_nton)[3:-2]:
            if comp_values[column] is not None:
                old_value = sum_comp_values.get(column, 0)
                sum_comp_values[column] = old_value + comp_values[column]
        #Append value for DoS to list, so we can calculate the mean afterwards
        if comp_values['direction of selection'] is not None:
            dos_list_a.append(comp_values['direction of selection'])

    #Calculate DoS average
    sum_comp_values['direction of selection'] = sum(dos_list_a) / len(dos_list_a)

    #Neutrality Index = Sum(X = Ds*Pn/(Ps+Ds)) / Sum(Y = Dn*Ps/(Ps+Ds))
    sum_comp_values['neutrality index'] = sum_comp_values['Ds*Pn/(Ps+Ds)'] / sum_comp_values['Dn*Ps/(Ps+Ds)']

    _append_statistics(calculations_file, 'total', sum_comp_values, sfs_max_nton)

def calculate_pnps(genome_ids_a, genome_ids_b, sico_files):
    """Calculate pN, pS, pN/pS & SFS values for all sico_files per clade, and write out statistics per SICO."""
    run_dir = tempfile.mkdtemp(prefix = 'calculate_pnps')

    #For each alignment create separate alignments for clade A & clade B genomes
    calculations_a_file = tempfile.mkstemp(suffix = '.tsv', prefix = 'calculations_a_')[1]
    calculations_b_file = tempfile.mkstemp(suffix = '.tsv', prefix = 'calculations_b_')[1]

    sfs_max_nton_a = len(genome_ids_a) / 2
    sfs_max_nton_b = len(genome_ids_b) / 2
    if 1 < len(genome_ids_a):
        _write_output_file_header(calculations_a_file, sfs_max_nton_a)
    if 1 < len(genome_ids_b):
        _write_output_file_header(calculations_b_file, sfs_max_nton_b)

    #Append all computed values to these lists, so we can compute sums and mean afterwards
    comp_values_a = []
    comp_values_b = []

    #Perform calculation for each sico_file
    for sico_file in sico_files:
        ali = AlignIO.read(sico_file, 'fasta')
        alignment_a = MultipleSeqAlignment(seqr for seqr in ali if seqr.id.split('|')[0] in genome_ids_a)
        alignment_b = MultipleSeqAlignment(seqr for seqr in ali if seqr.id.split('|')[0] in genome_ids_b)

        #1. gene name
        basename = os.path.split(sico_file)[1].split('.')[0]

        #Run codeml to calculate values for dn & ds
        codeml_file = run_codeml(create_directory(basename, inside_dir = run_dir), alignment_a, alignment_b)
        value_dict = parse_codeml_output(codeml_file)

        #Perform calculations for subaligments of each clade, if clade has more than one sequence; skipping outliers
        if 1 < len(alignment_a):
            comp_values = _perform_calculations(alignment_a, value_dict)
            comp_values_a.append(comp_values)
            _append_statistics(calculations_a_file, basename, comp_values, sfs_max_nton_a)
        if 1 < len(alignment_b):
            comp_values = _perform_calculations(alignment_b, value_dict)
            comp_values_b.append(comp_values)
            _append_statistics(calculations_b_file, basename, comp_values, sfs_max_nton_b)

    #"Finally, we might need a line which gives the sum for columns 3 to 15 plus the mean of column 16"
    if 1 < len(alignment_a):
        _append_sums_and_dos_average(calculations_a_file, sfs_max_nton_a, comp_values_a)
    if 1 < len(alignment_b):
        _append_sums_and_dos_average(calculations_b_file, sfs_max_nton_b, comp_values_b)

    #TODO Alternate method of calculating number of substitutions for independent X-axis of eventual graph

    #Clean up
    shutil.rmtree(run_dir)

    return calculations_a_file, calculations_b_file

#Using the standard NCBI Bacterial, Archaeal and Plant Plastid Code translation table (11).
BACTERIAL_CODON_TABLE = CodonTable.unambiguous_dna_by_id.get(11)

def _four_fold_degenerate_patterns():
    """Find patterns of 4-fold degenerate codons, wherein all third site substitutions code for the same amino acid."""
    letters = BACTERIAL_CODON_TABLE.nucleotide_alphabet.letters
    #Any combination of letters of length two
    for site12 in [''.join(prod) for prod in product(letters, repeat = 2)]:
        #4-fold when the length of the unique encoded amino acids for all possible third site nucleotides is exactly 1 
        if 1 == len(set([BACTERIAL_CODON_TABLE.forward_table.get(site12 + site3) for site3 in letters])):
            #Add regular expression pattern to the set of patterns
            yield '{0}[{1}]'.format(site12, letters)

FOUR_FOLD_DEGENERATE_PATTERNS = set(_four_fold_degenerate_patterns())

def _get_nton_name(nton, prefix = ''):
    """Given the number of strains in which a polymorphism/substitution is found, give the appropriate SFS name."""
    return prefix + {1:'singletons', 2:'doubletons', 3:'tripletons', 4:'quadrupletons', 5:'quintupletons'}.get(nton,
                                                                                                    str(nton) + '-tons')

def _perform_calculations(alignment, codeml_values):
    """Perform actual calculations on the alignment to determine pN, pS, SFS & the number of ignored cases per SICO."""
    synonymous_sfs = {}
    four_fold_syn_sfs = {}
    non_synonymous_sfs = {}
    four_fold_synonymous_sites = 0
    mixed_synonymous_polymorphisms = 0
    multiple_site_polymorphisms = 0

    #Calculate sequence_lengths here so we can handle alignments that are not multiples of three
    sequence_lengths = len(alignment[0]) - len(alignment[0]) % 3
    #Split into codon_alignments
    codon_alignments = [alignment[:, index:index + 3] for index in range(0, sequence_lengths, 3)]
    for codon_alignment in codon_alignments:
        #Get string representations of codons for simplicity 
        codons = [str(seqr.seq) for seqr in codon_alignment]

        #As per AEW: ignore codons with gaps, and codons with unresolved bases: Basically anything but ACGT
        if 0 < len(''.join(codons).translate(None, 'ACGTactg')):
            continue

        #Stop codons should have already been removed, but lets be sure anyway
        codons_nonstop = [codon for codon in codons if codon not in BACTERIAL_CODON_TABLE.stop_codons]

        #Retrieve translations of codons now that inconclusive & stop-codons have been removed
        translations = [BACTERIAL_CODON_TABLE.forward_table.get(codon) for codon in codons_nonstop]

        #Count unique translations across strains
        translation_usage = dict((aa, translations.count(aa)) for aa in set(translations))

        #Mutations are synonymous when all codons encode the same AA, and there are no skipped codons
        synonymous = len(translation_usage) == 1 and len(translations) == len(codon_alignment)

        #Retrieve nucleotides per site within the codon 
        site1 = [nucl for nucl in codon_alignment[:, 0]]
        site2 = [nucl for nucl in codon_alignment[:, 1]]
        site3 = [nucl for nucl in codon_alignment[:, 2]]

        #Count occurrences of distinct nucleotides across strains
        site1_usage = dict((nucl, site1.count(nucl)) for nucl in set(site1))
        site2_usage = dict((nucl, site2.count(nucl)) for nucl in set(site2))
        site3_usage = dict((nucl, site3.count(nucl)) for nucl in set(site3))

        #Sites are polymorphic if they contain more than one nucleotide
        site1_polymorphic = 1 < len(site1_usage)
        site2_polymorphic = 1 < len(site2_usage)
        site3_polymorphic = 1 < len(site3_usage)
        polymorphisms = site1_polymorphic, site2_polymorphic, site3_polymorphic

        #Continue with next codon if none of the sites is polymorphic
        if not any(polymorphisms):
            #But do increase the number of 4-fold synonymous sites if the pattern matches
            codon = codons_nonstop[0]
            for pattern in FOUR_FOLD_DEGENERATE_PATTERNS:
                if re.match(pattern, codon):
                    #Increase by one, as this site is for fold degenerate, even if it is not polymorphic
                    four_fold_synonymous_sites += 1
            continue

        #Determine if only one site is polymorphic by using boolean xor and not all
        single_site_polymorphism = site1_polymorphic ^ site2_polymorphic ^ site3_polymorphic and not all(polymorphisms)

        #Skip multiple site polymorphisms, but do keep a count of how many we encounter
        if not single_site_polymorphism:
            multiple_site_polymorphisms += 1
            continue

        #Determine which site_usage is the single site polymorphism
        polymorph_site_usage = site1_usage if site1_polymorphic else site2_usage if site2_polymorphic else site3_usage

        #Find the 'reference' nucleotide as (one of) the most occurring occupations in this site, so we can -1 later
        psu_values = polymorph_site_usage.values()
        reference_allele_count = max(psu_values)

        #Calculate the local site frequency spectrum, to be added to the gene-wide SFS later
        #We'll be using Site Frequency Spectrum to calculate the number of synonymous and non synonymous polymorphisms
        #Note: this requires a complete SFS across synonymous & non_synonymous polymorphisms, be careful when updating
        local_sfs = dict((ntimes, psu_values.count(ntimes)) for ntimes in set(psu_values))
        #Deduct one for the reference_allele_count, which should not count towards the SFS
        local_sfs[reference_allele_count] = local_sfs[reference_allele_count] - 1
        #Remove empty value as possible result of the above decrement operation
        if local_sfs[reference_allele_count] == 0:
            del local_sfs[reference_allele_count]

        def _update_sfs_with_local_sfs(sfs, local_sfs):
            """Add values from local_sfs to gene-wide sfs"""
            for maf, count in local_sfs.iteritems():
                prev_occupations = sfs.get(maf, 0)
                sfs[maf] = prev_occupations + count

        if synonymous:
            #If all polymorphisms encode for the same AA, we have multiple synonymous polymorphisms, where:
            #2 nucleotides = 1 polymorphism, 3 nucleotides = 2 polymorphisms, 4 nucleotides = 3 polymorphisms

            #Update synonymous SFS by adding values from local SFS
            _update_sfs_with_local_sfs(synonymous_sfs, local_sfs)

            #Codon is four fold degenerate if it matches any pattern in FOUR_FOLD_DEGENERATE_PATTERNS
            if site3_polymorphic:
                codon = codons_nonstop[0]
                for pattern in FOUR_FOLD_DEGENERATE_PATTERNS:
                    if re.match(pattern, codon):
                        #Update four fold degenerate SFS by adding values from local SFS
                        _update_sfs_with_local_sfs(four_fold_syn_sfs, local_sfs)
                        #Increase the number of four_fold synonymous sites here as well
                        four_fold_synonymous_sites += 1
        else: #not synonymous
            if len(polymorph_site_usage) == len(translation_usage):
                #If all polymorphisms encode for different AA, we have multiple non-synonymous polymorphisms, where:
                #2 nucleotides = 1 polymorphism, 3 nucleotides = 2 polymorphisms, 4 nucleotides = 3 polymorphisms

                #Update non synonymous SFS by adding values from local SFS
                _update_sfs_with_local_sfs(non_synonymous_sfs, local_sfs)
            else:
                #Some, but not all polymorphisms encode for different AA, making it unclear how this should be scored
                mixed_synonymous_polymorphisms += 1

    #Compute combined values from the above counted statistics 
    computed_values = _compute_values_from_statistics(len(alignment), sequence_lengths, codeml_values,
                                                      synonymous_sfs, non_synonymous_sfs, four_fold_syn_sfs,
                                                      four_fold_synonymous_sites)

    #Miscellaneous additional values
    computed_values['codons'] = sequence_lengths / 3
    computed_values['multiple site polymorphisms'] = multiple_site_polymorphisms
    computed_values['synonymous and non-synonymous polymorphisms mixed'] = mixed_synonymous_polymorphisms
    #Add COGs to output file
    computed_values['cogs'] = ','.join(find_cogs_in_sequence_records(alignment))

    return computed_values

def _compute_values_from_statistics(nr_of_strains, sequence_lengths, codeml_values,
                                    synonymous_sfs, non_synonymous_sfs, four_fold_syn_sfs, four_fold_synonymous_sites):
    """Compute values (mostly from the site frequency spectra) that we'll output according to provided formula's."""
    #12. number of non-synonymous sites for divergence from PAML
    #(this might be different to 3, because this will come from two randomly chosen sequences) = LnD
    paml_non_synonymous_sites = codeml_values['N']

    #14. number of synonymous sites from PAML = LsD
    paml_synonymous_sites = codeml_values['S']

    #Dictionary to hold computed values
    calc_values = {}

    #Add all values in codeml_values to calc_values
    calc_values.update(codeml_values)

    #2. number of strains (this will typically be the same for all genes)
    calc_values['strains'] = nr_of_strains

    def _add_sfs_values_in_columns(sfs_name, sfs):
        """Add values contained within SFS to named columns for singletons, doubletons, tripletons, etc.."""
        max_nton = nr_of_strains / 2 + 1
        for nton in range(1, max_nton):
            name = _get_nton_name(nton, prefix = sfs_name)
            value = sfs[nton] if nton in sfs else 0
            calc_values[name] = value

    #3. number of non-synonymous sites (see below for calculation) = LnP = L * LnD/(LnD+LsD)
    #where L is the length of the sequence from which the polymorphism data is taken (i.e. 3* no of codons used)
    calc_values['non-synonymous sites'] = sequence_lengths * paml_non_synonymous_sites / (paml_non_synonymous_sites +
                                                                                          paml_synonymous_sites)

    #4. SFS for non-synonymous polymorphsims
    _add_sfs_values_in_columns('non-synonymous sfs ', non_synonymous_sfs)

    #5. The sum of SFS for non-synonymous polymorphisms = Pn
    calc_values['non-synonymous polymorphisms'] = non_synonymous_polymorphisms = sum(non_synonymous_sfs.values())

    #6. number of synonymous sites (see below for calculation) = LsP = L * LsD/(LnD+LsD)
    #where L is the length of the sequence from which the polymorphism data is taken (i.e. 3* no of codons used)
    calc_values['synonymous sites'] = sequence_lengths * paml_synonymous_sites / (paml_non_synonymous_sites +
                                                                                  paml_synonymous_sites)

    #7. SFS for synonymous polymorphisms
    _add_sfs_values_in_columns('synonymous sfs ', synonymous_sfs)

    #8. The sum of the SFS for synonymous polymorphisms = Ps
    calc_values['synonymous polymorphisms'] = synonymous_polymorphisms = sum(synonymous_sfs.values())

    #9. number of 4-fold synonymous sites = L4  (all degenerate sites, _including_ those that are not polymorphic)
    calc_values['4-fold synonymous sites'] = four_fold_synonymous_sites

    #10. SFS for 4-fold synonymous polymorphisms
    _add_sfs_values_in_columns('4-fold synonymous sfs ', four_fold_syn_sfs)

    #11. Sum of the SFS for 4-folds = P4
    calc_values['4-fold synonymous polymorphisms'] = sum(four_fold_syn_sfs.values())

    #13. number of non-synonymous substitutions from PAML = Dn
    paml_non_synonymous_substitut = codeml_values['Dn']

    #15. number of synonymous substitutions from PAML = Ds
    paml_synonymous_substitutions = codeml_values['Ds']

    #16. Direction of selection = Dn/(Dn+Ds) - Pn/(Pn+Ps)
    paml_total_substitutions = paml_non_synonymous_substitut + paml_synonymous_substitutions
    total_polymorphisms = non_synonymous_polymorphisms + synonymous_polymorphisms
    #Prevent divide by zero by checking both values above are not null
    if paml_total_substitutions != 0 and total_polymorphisms != 0:
        calc_values['direction of selection'] = (paml_non_synonymous_substitut / paml_total_substitutions
                                                 - non_synonymous_polymorphisms / total_polymorphisms)
    else:
        #"the direction of selection is undefined if either Dn+Ds or Pn+Ps are zero":Not a number
        calc_values['direction of selection'] = None

    #These values will end up contributing to the Neutrality Index through NI = Sum(X) / Sum(Y)
    ps_plus_ds = (synonymous_polymorphisms + paml_synonymous_substitutions)
    if ps_plus_ds:
        #X = Ds*Pn/(Ps+Ds)
        calc_values['Ds*Pn/(Ps+Ds)'] = paml_synonymous_substitutions * non_synonymous_polymorphisms / ps_plus_ds
        #Y = Dn*Ps/(Ps+Ds)
        calc_values['Dn*Ps/(Ps+Ds)'] = paml_non_synonymous_substitut * synonymous_polymorphisms / ps_plus_ds
    else:
        calc_values['Ds*Pn/(Ps+Ds)'] = None
        calc_values['Dn*Ps/(Ps+Ds)'] = None

    return calc_values

def _get_column_headers_in_sequence(max_nton):
    """Return the column headers to the generated statistics files in the correct order, using max_nton for SFS max."""
    #Some initial values
    headers = ['cogs', 'strains', 'codons']

    def _write_sfs_column_names(prefix):
        """Return named columns for SFS singleton, doubleton, tripleton, etc.. upto max_nton."""
        for number in range(1, max_nton + 1):
            yield _get_nton_name(number, prefix)

    #Non synonymous
    headers.extend(['non-synonymous sites', 'non-synonymous polymorphisms'])
    headers.extend(_write_sfs_column_names('non-synonymous sfs '))
    #Synonymous
    headers.extend(['synonymous sites', 'synonymous polymorphisms'])
    headers.extend(_write_sfs_column_names('synonymous sfs '))
    #4-fold synonymous
    headers.extend(['4-fold synonymous sites', '4-fold synonymous polymorphisms'])
    headers.extend(_write_sfs_column_names('4-fold synonymous sfs '))
    #Miscellaneous additional statistics
    headers.append('multiple site polymorphisms')
    headers.append('synonymous and non-synonymous polymorphisms mixed')
    #PAML
    headers.extend(['N', 'Dn', 'S', 'Ds'])
    headers.append('Ds*Pn/(Ps+Ds)')
    headers.append('Dn*Ps/(Ps+Ds)')

    #Final values here are ignored when calculation sums over complete table
    headers.append('neutrality index')
    headers.append('direction of selection')

    return headers

def _write_output_file_header(calculations_file, max_nton):
    """Write header line for combined statistics file."""
    with open(calculations_file, mode = 'a') as append_handle:
        append_handle.write('ortholog\t')
        append_handle.write('\t'.join(_get_column_headers_in_sequence(max_nton)))
        append_handle.write('\n')

def _append_statistics(calculations_file, orthologname, comp_values, max_nton):
    """Append statistics for individual ortholog to genome-wide files."""
    with open(calculations_file, mode = 'a') as append_handle:
        append_handle.write(orthologname + '\t')
        append_handle.write('\t'.join(str(comp_values.get(column, ''))
                                      for column in _get_column_headers_in_sequence(max_nton)))
        append_handle.write('\n')

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: run_codeml.py 
--genomes-a=FILE    file with GenBank Project IDs from complete genomes table on each line for taxon A
--genomes-b=FILE    file with GenBank Project IDs from complete genomes table on each line for taxon B
--sico-zip=FILE     archive of aligned & trimmed single copy orthologous (SICO) genes
--table-a=FILE      destination file path for summary statistics table based on orthologs in taxon A
--table-b=FILE      destination file path for summary statistics table based on orthologs in taxon B
"""
    options = ['genomes-a', 'genomes-b', 'sico-zip', 'table-a', 'table-b']
    genome_a_ids_file, genome_b_ids_file, sico_zip, table_a, table_b = parse_options(usage, options, args)

    #Parse file containing GenBank GenBank Project IDs to extract GenBank Project IDs
    with open(genome_a_ids_file) as read_handle:
        genome_ids_a = [line.split()[0] for line in read_handle]
    with open(genome_b_ids_file) as read_handle:
        genome_ids_b = [line.split()[0] for line in read_handle]

    #Create run_dir to hold files relating to this run
    run_dir = tempfile.mkdtemp(prefix = 'run_codeml_')

    #Extract files from zip archive
    sico_files = extract_archive_of_files(sico_zip, create_directory('sicos', inside_dir = run_dir))

    #Actually do calculations
    tmp_table_tuple = calculate_pnps(genome_ids_a, genome_ids_b, sico_files)

    #Write the produced files to command line argument filenames
    shutil.move(tmp_table_tuple[0], table_a)
    shutil.move(tmp_table_tuple[1], table_b)

    #Remove unused files to free disk space 
    shutil.rmtree(run_dir)

    #Exit after a comforting log message
    log.info("Produced: \n%s\n%s", table_a, table_b)

if __name__ == '__main__':
    main(sys.argv[1:])

