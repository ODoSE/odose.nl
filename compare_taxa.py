#!/usr/bin/env python
"""Module to compare filtered and unfiltered taxa deduced from their respective trees to assert that they match."""

from divergence import parse_options
import logging as log
import sys

def fail(unfiltered_a, unfiltered_b, filtered_a, filtered_b):
    """Report error back to the user and exit with error code 1."""

    def _log_error_and_print_stderr(msg, dictionary = None):
        """Both log an error and print it to sys.stderr"""
        log.error(msg)
        print >> sys.stderr, msg
        if dictionary:
            for key, value in dictionary.iteritems():
                log.error('{0}\t{1}'.format(key, value))
                print >> sys.stderr, '{0}\t{1}'.format(key, value)

    _log_error_and_print_stderr('Unfiltered & filtered tree clusterings do not match!')
    _log_error_and_print_stderr('Unfiltered taxon A:', unfiltered_a)
    _log_error_and_print_stderr('Unfiltered taxon B:', unfiltered_b)
    _log_error_and_print_stderr('Filtered taxon A:', filtered_a)
    _log_error_and_print_stderr('Filtered taxon B:', filtered_b)

    sys.exit(1)

def main(args):
    """Main function called when run from command line or as part of pipeline."""
    usage = """
Usage: compare_taxa.py
--unfiltered-taxon-a=FILE    genome IDs for taxon A as deduced from phylogenetic tree of unfiltered concatemers
--unfiltered-taxon-b=FILE    genome IDs for taxon B as deduced from phylogenetic tree of unfiltered concatemers
--filtered-taxon-a=FILE      genome IDs for taxon A as deduced from phylogenetic tree of filtered concatemers
--filtered-taxon-b=FILE      genome IDs for taxon B as deduced from phylogenetic tree of filtered concatemers
"""
    options = ['unfiltered-taxon-a', 'unfiltered-taxon-b', 'filtered-taxon-a', 'filtered-taxon-b']
    unfiltered_a_file, unfiltered_b_file, filtered_a_file, filtered_b_file = parse_options(usage, options, args)

    #Parse file containing RefSeq project IDs to extract RefSeq project IDs
    with open(unfiltered_a_file) as read_handle:
        unfiltered_a = dict((line.split('\t')[0], line.strip().split('\t')[1]) for line in read_handle)
    with open(unfiltered_b_file) as read_handle:
        unfiltered_b = dict((line.split('\t')[0], line.strip().split('\t')[1]) for line in read_handle)
    with open(filtered_a_file) as read_handle:
        filtered_a = dict((line.split('\t')[0], line.strip().split('\t')[1]) for line in read_handle)
    with open(filtered_b_file) as read_handle:
        filtered_b = dict((line.split('\t')[0], line.strip().split('\t')[1]) for line in read_handle)

    #Otherwise fail after 
    if unfiltered_a.keys()[0] in filtered_a:
        if not (set(unfiltered_a.keys()) == set(filtered_a.keys())
                and set(unfiltered_b.keys()) == set(filtered_b.keys())):
            fail(unfiltered_a, unfiltered_b, filtered_a, filtered_b)
    else:
        if not (set(unfiltered_a.keys()) == set(filtered_b.keys())
                and set(unfiltered_b.keys()) == set(filtered_a.keys())):
            fail(unfiltered_a, unfiltered_b, filtered_b, filtered_a)

if __name__ == '__main__':
    main(sys.argv[1:])
