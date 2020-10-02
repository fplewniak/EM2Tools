"""
Some utilities for NGS data manipulation built on pysam and deeptools.
"""
#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
import multiprocessing
import pysam
import numpy


def mean_quality(ref, start, end, bam):
    """
    Computes the mean mapping quality at each position in the region of reference ref specified by start and end.

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :return: a list of mean mapping quality values along the region in the reference sequence
    """
    p = [0] * (end - start + 1)
    for col in bam.pileup(ref, start, end):
        if start <= col.reference_pos <= end:
            p[col.reference_pos - start] = numpy.mean(col.get_mapping_qualities())
    return p


def sum_quality(ref, start, end, bam):
    """
    Computes the sum of mapping quality values at each position in the region of reference ref specified by start
    and end.

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :return: a list of sums of mapping quality values along the region in the reference sequence
    """
    p = [0] * (end - start + 1)
    for col in bam.pileup(ref, start, end):
        if start <= col.reference_pos <= end:
            p[col.reference_pos - start] = numpy.sum(col.get_mapping_qualities())
    return p


def number_of_reads(ref, start, end, bam):
    """
    Counts the number of mapping reads at every position in the region of reference ref specified by start.

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :return: a list of number of reads mapping at every position along the region of the reference sequence
    """
    p = [0] * (end - start + 1)
    for col in bam.pileup(ref, start, end):
        if start <= col.reference_pos <= end:
            p[col.reference_pos - start] = col.nsegments
    return p


def position(ref, start, end, bam):
    """
    Returns the reference positions of every position in the specified ref region with at least on mapping read. This
    function may be useful only for testing purposes or display.

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :return: a list of positions in the reference sequence where there is at least a mapping read
    """
    p = [0] * (end - start + 1)
    for col in bam.pileup(ref, start, end):
        if start <= col.reference_pos <= end:
            p[col.reference_pos - start] = col.reference_pos
    return p


def mapping_qualities(ref, start, end, bam):
    """
    Returns all the mapping quality values in the specified reference region

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :return: a list of all the mapping quality values
    """
    qval = []
    for col in bam.pileup(ref, start, end):
        if start <= col.reference_pos <= end:
            qval += col.get_mapping_qualities()
    return qval


class MappingProfile:
    """
    Class for computing mapping quality or coverage profiles along reference sequences from mapped reads in a BAM file.
    Some methods are provided to compute different profiles but other functions may be used as well for specific
    purposes.
    """

    def __init__(self, bamfile, threads=multiprocessing.cpu_count()):
        self.bam = pysam.AlignmentFile(bamfile, threads=threads)

    def get_profiles(self, profile=mean_quality, references=None, start=0, end=None):
        """
        Method to get the profiles for all the requested references as a dictionary with reference names as keys and
        profiles as values.

        :param profile: the function responsible for computing profile values, should take a reference name, start and
         end positions, and the bam file handler, must return a profile
        :param references: a reference name or a list thereof
        :param start: start position for each reference sequence
        :param end: end position for each reference, if None, the whole reference sequence is used
        :return: the requested profile
        """
        if references is None:
            references = list(self.bam.references)
        if not isinstance(references, list):
            references = [references]
        profiles = {}
        for ref in references:
            ref_length = self.bam.get_reference_length(ref) - 1
            if start > ref_length:
                profiles[ref] = []
            else:
                stop = ref_length if end is None else min(end, ref_length)
                profiles[ref] = profile(ref, start, stop, self.bam)
        return profiles

    def get_values(self, profile=mapping_qualities, references=None, start=0, end=None):
        """
        Method to get all values of some mapping statistics in all references and regions pooled into one single array.

        :param profile: the function responsible for computing profile values, should take a reference name, start and
         end positions, and the bam file handler, must return a profile
        :param references: a reference name or a list thereof
        :param start: start position for each reference sequence
        :param end: end position for each reference, if None, the whole reference sequence is used
        :return: an array with all values
        """
        profiles = self.get_profiles(profile=profile, references=references, start=start, end=end)
        mapping_values = []
        for ref in profiles:
            mapping_values += profiles[ref]
        return mapping_values

    def get_references_with_max(self, profile=number_of_reads, references=None, start=0, end=None):
        """
        Returns the references whose profile contains the maximum value across profiles of all references

        :param profile: the function responsible for computing profile values, should take a reference name, start and
         end positions, and the bam file handler, must return a profile
        :param references: a reference name or a list thereof
        :param start: start position for each reference sequence
        :param end: end position for each reference, if None, the whole reference sequence is used
        :return: a tuple containing the maximum value and a list of all references reaching that value
        """
        profiles = self.get_profiles(profile=profile, references=references, start=start, end=end)
        max_value = 0
        ref_max = []
        for ref in profiles:
            if max(profiles[ref]) > max_value:
                max_value = max(profiles[ref])
                ref_max = [ref]
            elif max(profiles[ref]) == max_value:
                ref_max.append(ref)
        return max_value, ref_max
