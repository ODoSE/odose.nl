#!/usr/bin/env python
"""Module to generate a scatterplot of the ortholog trimming statistics."""

from rpy2 import robjects
from rpy2.robjects.vectors import DataFrame

__author__ = "Tim te Beek"
__contact__ = "brs@nbic.nl"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


def scatterplot(retained_threshold, trim_tuples):

    for trimmed_file, alignment_length, trimmed_length, percentage in trim_tuples:
        print trimmed_file, alignment_length, trimmed_length, percentage

    #Extract desired columns
    original_lengths = [trim_tuple[1] for trim_tuple in trim_tuples]
    percentages_retained = [trim_tuple[3] for trim_tuple in trim_tuples]
    row_retained = [trim_tuple[3] >= retained_threshold for trim_tuple in trim_tuples]

    #Convert to R objects
    original_lengths = robjects.IntVector(original_lengths)
    percentages_retained = robjects.FloatVector(percentages_retained)
    row_retained = robjects.BoolVector(row_retained)


    dataframe = DataFrame({
                           'original length': original_lengths,
                           'percentage retained': percentages_retained,
                           #'row retained': row_retained
                           })
    print dataframe

    #Try ggplot2
    #from rpy2.robjects.lib import ggplot2
    #gp = ggplot2.ggplot(dataframe)

    #pp = gp + \
    #     ggplot2.aes_string(x='original length', y='percentage retained') + \
    #     ggplot2.geom_point()
    #
    #pp.plot()


    #Wait before exiting
    #raw_input()


