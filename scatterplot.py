#!/usr/bin/env python
"""Module to generate a scatterplot of the ortholog trimming statistics."""

from rpy2 import robjects

__author__ = "Tim te Beek"
__copyright__ = "Copyright 2011, Netherlands Bioinformatics Centre"
__license__ = "MIT"


def scatterplot(retained_threshold, trim_tuples, scatterplot_file):
    '''
    Create a scatterplot for the given retained percentage and statistics using R & Rpy2 on provided path.
    @param retained_threshold: the percentage below which sequences are filtered 
    @param trim_tuples: trim statistics as written to or read from trim stats file
    @param scatterplot_file: destination output path for created scatterplot PDF
    '''
    # Extract desired columns
    original_lengths = [trim_tuple[1] for trim_tuple in trim_tuples]
    percentages_retained = [trim_tuple[3] for trim_tuple in trim_tuples]
    row_retained = ['green' if trim_tuple[3] >= retained_threshold
                    else 'red' if trim_tuple[3] > 0 else 'black'
                    for trim_tuple in trim_tuples]

    # Convert to R objects
    original_lengths = robjects.IntVector(original_lengths)
    percentages_retained = robjects.FloatVector(percentages_retained)
    row_retained = robjects.StrVector(row_retained)

    # Open pdf plotter
    from rpy2.robjects.packages import importr
    grdevices = importr('grDevices')
    grdevices.pdf(file=scatterplot_file)  # pylint: disable=E1101

    # Plot using standard plot
    r = robjects.r
    r.plot(x=percentages_retained,
           y=original_lengths,
           col=row_retained,
           main='Sequences retained in alignment & trimming',
           xlab='Percentage retained',
           ylab='Original length',
           pch=18)

    # Add legend
    r.legend('top',
             legend=['retained', 'filtered', 'too large indel'],
             fill=robjects.StrVector(['green', 'red', 'black']))

    # Close the graphical device
    grdevices.dev_off()  # pylint: disable=E1101
