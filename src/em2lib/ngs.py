"""
Some utilities for NGS data manipulation built on pysam and deeptools.
"""
#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
import multiprocessing
import warnings

from pysam import AlignmentFile
import numpy
import gffpandas.gffpandas as gffpd
from pandas import DataFrame
from em2lib.table import get_max
from em2lib.table import select_rows
from em2lib.gff import select_features


def mean_quality(ref, start, end, bam, **kwargs):
    """
    Computes the mean mapping quality at each position in the region of reference ref specified by start and end.

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :param **kwargs: supplementary keyword parameters that will be passed to the pileup() method, allowing to filter
     reads by flag, quality, etc.
    :return: a list of mean mapping quality values along the region in the reference sequence
    """
    p = [0] * (end - start + 1)
    for col in bam.pileup(ref, start, end, **kwargs):
        if start <= col.reference_pos <= end:
            p[col.reference_pos - start] = numpy.mean(col.get_mapping_qualities())
    return p


def sum_quality(ref, start, end, bam, **kwargs):
    """
    Computes the sum of mapping quality values at each position in the region of reference ref specified by start
    and end.

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :param **kwargs: supplementary keyword parameters that will be passed to the pileup() method, allowing to filter
     reads by flag, quality, etc.
    :return: a list of sums of mapping quality values along the region in the reference sequence
    """
    p = [0] * (end - start + 1)
    for col in bam.pileup(ref, start, end, **kwargs):
        if start <= col.reference_pos <= end:
            p[col.reference_pos - start] = numpy.sum(col.get_mapping_qualities())
    return p


def number_of_reads(ref, start, end, bam, **kwargs):
    """
    Counts the number of mapping reads at every position in the region of reference ref specified by start.

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :param **kwargs: supplementary keyword parameters that will be passed to the pileup() method, allowing to filter
     reads by flag, quality, etc.
    :return: a list of number of reads mapping at every position along the region of the reference sequence
    """
    p = [0] * (end - start + 1)
    for col in bam.pileup(ref, start, end, **kwargs):
        if start <= col.reference_pos <= end:
            p[col.reference_pos - start] = col.nsegments
    return p


def position(ref, start, end, bam, **kwargs):
    """
    Returns the reference positions of every position in the specified ref region with at least on mapping read. This
    function may be useful only for testing purposes or display.

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :param **kwargs: supplementary keyword parameters that will be passed to the pileup() method, allowing to filter
     reads by flag, quality, etc.
    :return: a list of positions in the reference sequence where there is at least a mapping read
    """
    p = [0] * (end - start + 1)
    for col in bam.pileup(ref, start, end, **kwargs):
        if start <= col.reference_pos <= end:
            p[col.reference_pos - start] = col.reference_pos
    return p


def mapping_qualities(ref, start, end, bam, **kwargs):
    """
    Returns all the mapping quality values in the specified reference region

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :param **kwargs: supplementary keyword parameters that will be passed to the pileup() method, allowing to filter
     reads by flag, quality, etc.
    :return: a list of all the mapping quality values
    """
    qval = []
    for col in bam.pileup(ref, start, end, **kwargs):
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
        self.bam = AlignmentFile(bamfile, threads=threads)

    def get_profiles(self, profile=mean_quality, references=None, start=0, end=None, gff=None, ftype=None, **kwargs):
        """
        Method to get the profiles for all the requested references as a dictionary with reference names as keys and
        profiles as values.

        :param profile: the function responsible for computing profile values, should take a reference name, start and
         end positions, and the bam file handler, must return a profile
        :param references: a reference name or a list thereof for reference sequences to compute profiles of
        :param start: start position for each reference sequence (0-based). If a GFF file is passed, then this value
         is relative to the start of each feature.
        :param end: end position for each reference, if None, the whole reference sequence is used (0-based). If a GFF
         file is passed, then this value is relative to the start of each feature.
        :param gff: GFF file name defining features found in reference sequences
        :param ftype: the type of features to compute the statistics of. May be a list for multiple types
        :return: a dictionary of dictionaries containing the profile for each selected feature or reference sequence
         and the associated information (feature location, region which the profile was computed over, strand, etc.)
        """
        if gff is not None:
            # if a GFF file was registered, select features according to references and type
            gff_df = select_features(gff, references=references, ftype=ftype)
        else:
            # else, if no GFF file was registered, then take references sequences
            # ensure references is None or a list
            if references is not None and not isinstance(references, list):
                references = [references]
            gff_df = DataFrame([{'seq_id': ref, 'start': 1, 'end': self.bam.get_reference_length(ref), 'strand': '0',
                                 'locus_tag': ref}
                                for ref in self.bam.references if references is None or ref in references])
        profiles = {}
        # for each row in DataFrame containing annotations
        for index, row in gff_df.iterrows():
            # shift positions relatively to the feature location
            begin = row['start'] + start - 1
            stop = row['end'] - 1 if end is None else min(row['end'] - 1, row['start'] + end - 1)
            # just to make sure there is a locus_tag to include as a key
            tag = row['seq_id'] + ':' + str(begin) + ':' + str(stop) if row['locus_tag'] is None else row['locus_tag']
            # compute the requested profile along the reference sequence or feature or region thereof
            profiles[tag] = {'seq_id': row['seq_id'], 'start': row['start'], 'end': row['end'], 'strand': row['strand'],
                             'from': begin + 1, 'to': stop + 1}
            profiles[tag]['profile'] = profile(row['seq_id'], begin, stop, self.bam, **kwargs) if begin <= stop else []
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


class MappingStatistics:
    """
    A class to compute a variety of read mapping statistics from a bam file and optionally from annotations in a GFF
    file.
    """

    def __init__(self, bamfile, gff=None, threads=multiprocessing.cpu_count()):
        self.bam = AlignmentFile(bamfile, threads=threads)
        self.mapping_profile = MappingProfile(bamfile=bamfile, threads=threads)
        self.gff = gff

    def get_statistics(self, profile=mean_quality, references=None, start=0, end=None, ftype=None, **kwargs):
        """
        Return mapping statistics according to the profile function computed on the specified reference sequences,
        positions and type of features. If no GFF file has been passed to the current instance, then complete sequences
        are used. If a GFF file has been defined, then the positions apply relative to feature locations. This may be
        useful to get statistics for instance about the first 10 bases of every feature, etc.

        :param profile: the function responsible for computing profile values, should take a reference name, start and
         end positions, and the bam file handler, must return a profile
        :param references: a reference name or a list thereof
        :param start: start position for each reference sequence or feature. If a GFF file is passed, then this value
         is relative to the start of each feature.
        :param end: end position for each reference or feature, if None, this is equivalent to the end of the reference
         sequence or feature. If a GFF file is passed, then this value is relative to the start of each feature.
        :param ftype: the type of features to compute the statistics of.
        :return:
        """
        # prepare the DataFrame to hold the statistics
        stat_df = DataFrame(columns=['feature', 'ref_id', 'start', 'end', 'strand',
                                     'from', 'to', 'mean', 'median', 'min', 'max'])
        profiles = self.mapping_profile.get_profiles(profile=profile, references=references, start=start, end=end,
                                                     gff=self.gff, ftype=ftype, **kwargs)
        for feature in profiles:
            prof = profiles[feature]
            stat_df = stat_df.append(DataFrame([{'feature': feature, 'ref_id': prof['seq_id'],
                                                 'start': prof['start'], 'end': prof['end'], 'strand': prof['strand'],
                                                 'from': prof['from'], 'to': prof['to'],
                                                 'total': numpy.sum(prof['profile']),
                                                 'mean': numpy.mean(prof['profile']),
                                                 'median': numpy.median(prof['profile']),
                                                 'min': numpy.min(prof['profile']),
                                                 'max': numpy.max(prof['profile'])}]), ignore_index=True)
        return stat_df

    def select_refs(self, profile=number_of_reads, references=None, start=0, end=None, ftype=None, query=get_max,
                    **kwargs):
        """
        Returns sequences (references or features) whose statistics correspond to a given criterion. This method is
        actually a shortcut, using successively the get_statistics() method above and the table.select_rows() function.
        If several row selections on the same dataset are required, it would be more efficient to store the statistics
        in a DataFrame and then select the rows as many times as necessary in order to avoid recomputing the statistics
        every time.

        :param profile: the function responsible for computing profile values, should take a reference name, start and
         end positions, and the bam file handler, must return a profile
        :param references: a reference name or a list thereof
        :param start: start position for each reference sequence or feature. If a GFF file is passed, then this value
         is relative to the start of each feature.
        :param end: end position for each reference or feature, if None, this is equivalent to the end of the reference
         sequence or feature. If a GFF file is passed, then this value is relative to the start of each feature.
        :param ftype: the type of features to compute the statistics of.
        :param query: the selection function (returning a boolean value for each row according to the implemented test)
         or a string representing a pandas.DataFrame query on columns such as 'max - min <= 10' to select rows where
         difference between max and min is lesser or equal than 10
        :param **kwargs: supplementary key arguments passed to the function specified in query if not a string.
        :return: a tuple containing the maximum value and a list of all references reaching that value
        """
        warnings.warn('The select_refs method is deprecated and will be removed. '
                      'Create a statistics DataFrame first then use table.select_rows() method on that DataFrame.',
                      DeprecationWarning, stacklevel=2)
        stat_df = self.get_statistics(profile=profile, references=references, start=start, end=end, ftype=ftype)
        return select_rows(stat_df, query=query, **kwargs)
