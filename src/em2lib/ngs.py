"""
Some utilities for NGS data manipulation built on pysam and deeptools.
"""
#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
import multiprocessing
from pysam import AlignmentFile
import numpy
import gffpandas.gffpandas as gffpd
from pandas import DataFrame
from em2lib.table import get_max
from em2lib.table import select_rows


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
        self.bam = AlignmentFile(bamfile, threads=threads)

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


class MappingStatistics:
    """
    A class to compute a variety of read mapping statistics from a bam file and optionally from annotations in a GFF
    file.
    """

    def __init__(self, bamfile, gff=None, threads=multiprocessing.cpu_count()):
        self.bam = AlignmentFile(bamfile, threads=threads)
        self.mapping_profile = MappingProfile(bamfile=bamfile, threads=threads)
        self.gff = gffpd.read_gff3(gff) if gff is not None else None

    def get_statistics(self, profile=mean_quality, references=None, start=0, end=None, ftype=None):
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
        # ensure type is None or a list (required by filter_feature_of_type(type))
        if ftype is not None and not isinstance(ftype, list):
            ftype = [ftype]
        # ensure references is None or a list
        if references is not None and not isinstance(references, list):
            references = [references]
        if self.gff is not None:
            # if a GFF file was registered, then take references sequences
            gff_df = self.gff.attributes_to_columns() if ftype is None \
                else self.gff.filter_feature_of_type(ftype).attributes_to_columns()
            if references is not None:
                gff_df = gff_df.loc[gff_df['seq_id'].isin(references)]
            gff_df = gff_df.astype({'start': int, 'end': int})
        else:
            # else, if no GFF file was registered, then take features as specified by ftype
            gff_df = DataFrame([{'seq_id': ref, 'start': 1, 'end': self.bam.get_reference_length(ref), 'locus_tag': ref}
                                for ref in self.bam.references if references is None or ref in references])
        # prepare the DataFrame to hold the statistics
        stat_df = DataFrame(columns=['feature', 'ref_id', 'start', 'end', 'from', 'to', 'mean', 'median', 'min', 'max'])
        # for each row in DataFrame containing annotations
        for index, row in gff_df.iterrows():
            # shift positions relatively to the feature location
            begin = row['start'] + start - 1
            stop = row['end'] - 1 if end is None else min(row['end'], row['start'] + end - 1)
            # compute the requested values along the reference sequence or feature or region thereof
            prof = self.mapping_profile.get_profiles(profile=profile, references=row['seq_id'], start=begin, end=stop)
            # just to make sure there is a locus_tag to include in the DataFrame
            if 'locus_tag' not in row:
                row['locus_tag'] = row['seq_id']
            # fill the DataFrame with the values for the current feature, start and end are the feature's coordinates
            # while from and to are the coordinates of the region used to compute the statistics
            stat_df = stat_df.append(DataFrame([{'feature': row['locus_tag'], 'ref_id': row['seq_id'],
                                                 'start': row['start'], 'end': row['end'],
                                                 'from': begin + 1, 'to': stop + 1,
                                                 'mean': numpy.mean(prof[row['seq_id']]),
                                                 'median': numpy.median(prof[row['seq_id']]),
                                                 'min': numpy.min(prof[row['seq_id']]),
                                                 'max': numpy.max(prof[row['seq_id']])}]), ignore_index=True)
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
        :param **kwargs: suplementary key arguments passed to the function specified in query if not a string.
        :return: a tuple containing the maximum value and a list of all references reaching that value
        """
        stat_df = self.get_statistics(profile=profile, references=references, start=start, end=end, ftype=ftype)
        return select_rows(stat_df, query=query, **kwargs)
