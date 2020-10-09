"""
Some utilities for NGS data manipulation built on pysam and deeptools.
"""
#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
import multiprocessing
import warnings

import numpy
from pandas import DataFrame
from pysam import AlignmentFile

from em2lib.gff import select_features
from em2lib.table import get_max
from em2lib.table import select_rows


def read_is_first_forward(read):
    """
    check_read(read) type of function to check read is first in pair and maps on the forward strand

    :param read: the read to check
    :return: True if read is first in pair and maps on the forward strand
    """
    return read.is_read1 and not read.is_reverse


def read_is_first_reverse(read):
    """
    check_read(read) type of function to check read is first in pair and maps on the reverse strand

    :param read: the read to check
    :return: True if read is first in pair and maps on the reverse strand
    """
    return read.is_read1 and read.is_reverse


def read_is_second_forward(read):
    """
    check_read(read) type of function to check read is second in pair and maps on the forward strand

    :param read: the read to check
    :return: True if read is second in pair and maps on the forward strand
    """
    return read.is_read2 and not read.is_reverse


def read_is_second_reverse(read):
    """
    check_read(read) type of function to check read is second in pair and maps on the reverse strand

    :param read: the read to check
    :return: True if read is second in pair and maps on the reverse strand
    """
    return read.is_read2 and read.is_reverse


def read_is_unpaired_fwd(read):
    """
    check_read(read) type of function to check read is not paired and maps on the forward strand

    :param read: the read to check
    :return: True if read is not paired and maps on the forward strand
    """
    return not read.is_paired and not read.is_reverse


def read_is_unpaired_rev(read):
    """
    check_read(read) type of function to check read is not paired and maps on the reverse strand

    :param read: the read to check
    :return: True if read is not paired and maps on the reverse strand
    """
    return not read.is_paired and read.is_reverse


def mean_quality(ref, start, end, bam, **kwargs):
    """
    Computes the mean mapping quality at each position in the region of reference ref specified by start and end.

    kwargs: supplementary keyword parameters that will be passed to the pileup() method, allowing to filter
    reads by flag, quality, etc.

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :return: a list of mean mapping quality values along the region in the reference sequence
    """
    prof = [0] * (end - start + 1)
    for col in bam.pileup(ref, start, end, **kwargs):
        if start <= col.reference_pos <= end:
            prof[col.reference_pos - start] = numpy.mean(col.get_mapping_qualities())
    return prof


def sum_quality(ref, start, end, bam, **kwargs):
    """
    Computes the sum of mapping quality values at each position in the region of reference ref specified by start
    and end.

    kwargs: supplementary keyword parameters that will be passed to the pileup() method, allowing to filter
    reads by flag, quality, etc.

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :return: a list of sums of mapping quality values along the region in the reference sequence
    """
    prof = [0] * (end - start + 1)
    for col in bam.pileup(ref, start, end, **kwargs):
        if start <= col.reference_pos <= end:
            prof[col.reference_pos - start] = numpy.sum(col.get_mapping_qualities())
    return prof


def number_of_reads(ref, start, end, bam, **kwargs):
    """
    Counts the number of mapping reads at every position in the region of reference ref specified by start.

    kwargs: supplementary keyword parameters that will be passed to the pileup() method, allowing to filter
    reads by flag, quality, etc.

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :return: a list of number of reads mapping at every position along the region of the reference sequence
    """
    prof = [0] * (end - start + 1)
    for col in bam.pileup(ref, start, end, **kwargs):
        if start <= col.reference_pos <= end:
            prof[col.reference_pos - start] = col.nsegments
    return prof


def position(ref, start, end, bam, **kwargs):
    """
    Returns the reference positions of every position in the specified ref region with at least on mapping read. This
    function may be useful only for testing purposes or display.

    kwargs: supplementary keyword parameters that will be passed to the pileup() method, allowing to filter
    reads by flag, quality, etc.

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
    :return: a list of positions in the reference sequence where there is at least a mapping read
    """
    prof = [0] * (end - start + 1)
    for col in bam.pileup(ref, start, end, **kwargs):
        if start <= col.reference_pos <= end:
            prof[col.reference_pos - start] = col.reference_pos
    return prof


def mapping_qualities(ref, start, end, bam, **kwargs):
    """
    Returns all the mapping quality values in the specified reference region

    kwargs: supplementary keyword parameters that will be passed to the pileup() method, allowing to filter
    reads by flag, quality, etc.

    :param ref: reference sequence in bam file
    :param start: 0-based start of region to compute the profile of
    :param end: 0-based end of region to compute the profile of (this position is included in the profile)
    :param bam: pysam AlignmentFile handler
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

    def count_reads(self, ref, start, end, mapq=0, read_callback='nofilter'):
        """
        Counts reads in a reference sequence or a region. Reads can be filtered by mapping quality and using a
        read_callback function or instruction ('nofilter' or 'all')

        :param ref: the reference sequence
        :param start: start position in the reference sequence
        :param end: end position in the reference sequence
        :param mapq: the minimum mapping quality required to count a read
        :param read_callback: 'nofilter', 'all' or a function of the type check_read(read) returning True if a read must
         be counted. This parameter is passed to the pysam.AlignmentFile.count() method to filter reads.
        :return: the count of reads mapping on the reference sequence or the region
        """
        if read_callback == 'nofilter':
            def callback(x):
                """
                callback function returning True if mapping quality is greater than the threshold
                """
                return x.mapping_quality >= mapq
        elif read_callback == 'all':
            def callback(x):
                """
                callback function returning True if mapping quality is greater than the threshold and none of the
                following flags are set: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
                """
                check_all_callback = not (x.is_unmapped or x.is_duplicate or x.is_secondary or x.is_qcfail)
                return x.mapping_quality >= mapq and check_all_callback
        else:
            def callback(x):
                """
                callback function returning True if mapping quality is greater than the threshold and the specified
                callback function returns True
                """
                return read_callback(x) and x.mapping_quality >= mapq
        return self.bam.count(contig=ref, start=start, stop=end, read_callback=callback)

    def strand_statistics(self, references=None, start=0, end=None, ftype=None, mapq=0):
        """
        Counts first and second reads mapping on reverse and forward strands.

        :param references: a reference name or a list thereof
        :param start: start position for each reference sequence or feature. If a GFF file is passed, then this value
         is relative to the start of each feature.
        :param end: end position for each reference or feature, if None, this is equivalent to the end of the reference
         sequence or feature. If a GFF file is passed, then this value is relative to the start of each feature.
        :param ftype: the type of features to compute the statistics of.
        :param mapq: the minimum mapping quality value to count a read
        :return: a DataFrame containing the statistics
        """
        stat_df = DataFrame(columns=['feature', 'type', 'ref_id', 'start', 'end', 'strand',
                                     'from', 'to', '1st rev', '2nd fwd', '1st fwd', '2nd rev', 'total'])

        if self.gff is not None:
            # if a GFF file was registered, select features according to references and type
            gff_df = select_features(self.gff, references=references, ftype=ftype)
        else:
            # else, if no GFF file was registered, then take references sequences
            # ensure references is None or a list
            if references is not None and not isinstance(references, list):
                references = [references]
            gff_df = DataFrame([{'seq_id': ref, 'start': 1, 'end': self.bam.get_reference_length(ref), 'strand': '0',
                                 'locus_tag': ref, 'type': 'source'}
                                for ref in self.bam.references if references is None or ref in references])

        for index, row in gff_df.iterrows():
            # shift positions relatively to the feature location
            begin = row['start'] + start - 1
            stop = row['end'] - 1 if end is None else min(row['end'] - 1, row['start'] + end - 1)
            # just to make sure there is a locus_tag to include as a key
            tag = row['seq_id'] + ':' + str(begin) + ':' + str(stop) if row['locus_tag'] is None else row['locus_tag']
            fr = self.count_reads(row['seq_id'], begin, end, read_callback=read_is_first_reverse, mapq=mapq)
            ff = self.count_reads(row['seq_id'], begin, end, read_callback=read_is_first_forward, mapq=mapq)
            sr = self.count_reads(row['seq_id'], begin, end, read_callback=read_is_second_reverse, mapq=mapq)
            sf = self.count_reads(row['seq_id'], begin, end, read_callback=read_is_second_forward, mapq=mapq)
            total = fr + ff + sf + sr
            stat_df = stat_df.append(DataFrame([{'feature': tag, 'type': row['type'], 'ref_id': row['seq_id'],
                                                 'start': row['start'], 'end': row['end'], 'strand': row['strand'],
                                                 'from': begin + 1, 'to': stop + 1,
                                                 '1st rev': fr,
                                                 '2nd fwd': sf,
                                                 '1st fwd': ff,
                                                 '2nd rev': sr,
                                                 'total': total}]), ignore_index=True)
        return stat_df

    def get_statistics(self, profile=mean_quality, references=None, start=0, end=None, ftype=None, **kwargs):
        """
        Return mapping statistics according to the profile function computed on the specified reference sequences,
        positions and type of features. If no GFF file has been passed to the current instance, then complete sequences
        are used. If a GFF file has been defined, then the positions apply relative to feature locations. This may be
        useful to get statistics for instance about the first 10 bases of every feature, etc.

        kwargs are keyword arguments passed to the get_profiles() method

        :param profile: the function responsible for computing profile values, should take a reference name, start and
         end positions, and the bam file handler, must return a profile
        :param references: a reference name or a list thereof
        :param start: start position for each reference sequence or feature. If a GFF file is passed, then this value
         is relative to the start of each feature.
        :param end: end position for each reference or feature, if None, this is equivalent to the end of the reference
         sequence or feature. If a GFF file is passed, then this value is relative to the start of each feature.
        :param ftype: the type of features to compute the statistics of.
        :return: a DataFrame containing the statistics
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
        Deprecated, use get_statistics() method above then the table.select_rows() function instead.

        Returns sequences (references or features) whose statistics correspond to a given criterion. This method is
        actually a shortcut, using successively the get_statistics() method above and the table.select_rows() function.
        If several row selections on the same dataset are required, it would be more efficient to store the statistics
        in a DataFrame and then select the rows as many times as necessary in order to avoid recomputing the statistics
        every time.

        kwargs: supplementary key arguments passed to the function specified in query if not a string.

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
        :return: a tuple containing the maximum value and a list of all references reaching that value
        """
        warnings.warn('The select_refs method is deprecated and will be removed. '
                      'Create a statistics DataFrame first then use table.select_rows() method on that DataFrame.',
                      DeprecationWarning, stacklevel=2)
        stat_df = self.get_statistics(profile=profile, references=references, start=start, end=end, ftype=ftype)
        return select_rows(stat_df, query=query, **kwargs)
