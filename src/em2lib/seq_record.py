"""
Extension module to the Biopython Bio.SeqRecord module
"""

#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
import warnings

from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from em2lib.seq import SeqEM2
from em2lib.seq_feature import SeqFeatureEM2


class SeqRecordEM2(SeqRecord):
    """
    Extension to Biopython SeqRecord class
    """

    def __init__(self, seq, **kwargs):
        super().__init__(seq, **kwargs)

    def __lt__(self, other):
        return self.id.__lt__(other.id)

    def __le__(self, other):
        return self.id.__le__(other.id)

    def __le___(self, other):
        return self.__le__(other)

    def __eq__(self, other):
        return self.id.__eq__(other.id)

    def __ne__(self, other):
        return self.id.__ne__(other.id)

    def __gt__(self, other):
        return self.id.__gt__(other.id)

    def __ge__(self, other):
        return self.id.__ge__(other.id)

    def overlap(self, start, end=None, strand=0):
        """
        Retrieves features that overlap a given position range.

        :param strand: strand specification of features to be returned. If strand is 0, then
        features on both strands are returned. If feature.strand is 0, then all strands will match.
        :param start: start of range
        :param end: end of range, if None, then end=start

        :return: a list of overlaping features
        """
        if end is None:
            end = start
        return [f for f in self.features if f.overlaps(start, end) and f.strand * strand >= 0]

    def feature_after(self, position, strand=0, nearest=False):
        """
        Retrieves the features immediately after (but not overlaping) the specified position, on one
        strand or both.
        If nearest is True, then only the nearest ones are returned.

        :param nearest: if True, only the nearest features are returned. This only makes sense when
        strand is O
        :param strand: strand specification of features to be returned. If strand is 0, then
        features on both strands are returned
        :param position: the position

        :return: a list of features after the specified position
        """
        after_fwd = {f: f.location.start for f in self.features if
                     f.location.start > position and f.strand in [0, 1]}
        after_rev = {f: f.location.end for f in self.features if
                     f.location.end < position and f.strand in [0, -1]}
        fwd_list = [k for k, v in after_fwd.items() if v == min(after_fwd.values())]
        rev_list = [k for k, v in after_rev.items() if v == max(after_rev.values())]
        if strand == 1:
            return fwd_list
        if strand == -1:
            return rev_list
        if nearest is True:
            if len(after_fwd) == len(after_rev) == 0:
                nearest_list = []
            elif len(after_fwd) == 0:
                nearest_list = rev_list
            elif len(after_rev) == 0:
                nearest_list = fwd_list
            elif fwd_list[0].location.start - position < position - rev_list[0].location.end:
                return fwd_list
            else:
                nearest_list = rev_list
            return nearest_list
        return list(set(fwd_list + rev_list))

    def feature_before(self, position, strand=0, nearest=False):
        """
        Retrieves the features immediately before (but not overlaping) the specified position, on
        one strand or both.
        If nearest is True, then only the nearest ones are returned.

        :param nearest: if True, only the nearest features are returned. This only makes sense when
        strand is O
        :param strand: strand specification of features to be returned. If strand is 0, then
        features on both strands are returned
        :param position: the position

        :return: a list of features before the specified position
        """
        before_fwd = {f: f.location.end for f in self.features if
                      f.location.end < position and f.strand in [0, 1]}
        before_rev = {f: f.location.start for f in self.features if
                      f.location.start > position and f.strand in [0, -1]}
        fwd_list = [k for k, v in before_fwd.items() if v == max(before_fwd.values())]
        rev_list = [k for k, v in before_rev.items() if v == min(before_rev.values())]
        if strand == 1:
            return fwd_list
        if strand == -1:
            return rev_list
        if nearest is True:
            if len(before_fwd) == len(before_rev) == 0:
                nearest_list = []
            elif len(before_fwd) == 0:
                nearest_list = rev_list
            elif len(before_rev) == 0:
                nearest_list = fwd_list
            elif position - fwd_list[0].location.end < rev_list[0].location.start - position:
                nearest_list = fwd_list
            else:
                nearest_list = rev_list
            return nearest_list
        return list(set(fwd_list + rev_list))

    def surrounding_features(self, position, strand=0, nearest=False):
        """
        Retrieves all the features around a given position but not overlaping it. If nearest is
        True, then only the nearest features are returned.

        :param position: the position
        :param nearest: if True, only the nearest features are returned.
        :param strand: strand specification of features to be returned. If strand is 0, then
        features on both strands are returned

        :return: a list of features around the specified position
        """
        if nearest is False:
            return list(
                set(self.feature_after(position, strand) + self.feature_before(position, strand)))
        feat_list = {f: min(abs(position - f.location.start), abs(position - f.location.end)) for f
                     in list(set(self.feature_after(position, strand, True)
                                 + self.feature_before(position, strand, True)))}
        return [f for f, v in feat_list.items() if v == min(feat_list.values())]

    def join(self, other=None, offset=0, keepself=True):
        """
        Joins two SeqRecordEM2 objects into a new one representing the resulting merged sequence

        :param keepself: if True and overlapping subsequences are different, then keep sequence from
         self record, otherwise keep the sequence of other record.
        :param other: the other SeqRecordEM2 object
        :param offset: the offset of the two sequences. If the value is negative, then the two
        sequences overlap.

        :return: the result of merging records as a new SeqRecordEM2 object
        """
        if len(self.seq) + offset < 0:
            return other.join(self, offset=-len(self.seq) - len(other.seq) - offset,
                              keepself=not keepself)

        if offset >= 0:
            new_seq = str(self.seq) + self.seq.any_residue * offset + str(other.seq)
        else:
            if str(self.seq)[offset:] != str(other.seq)[0:-offset]:
                warnings.warn('Warning!!! Overlapping subsequences are different.')
            if keepself is True:
                new_seq = str(self.seq) + str(other.seq)[-offset:]
            else:
                new_seq = str(self.seq)[0:offset] + str(other.seq)

        if self.seq.is_protein() and other.seq.is_protein():
            new_record = SeqRecordEM2(SeqEM2.protein(new_seq))
        elif not (self.seq.is_protein() or other.seq.is_protein()):
            new_record = SeqRecordEM2(SeqEM2.dna(new_seq))
        else:
            raise ValueError('Sequences are not of the same type. It is impossible to join them.')

        for feature in self.features:
            new_record.features.append(
                SeqFeatureEM2(parent=new_record, location=feature.location, strand=feature.strand,
                              id=feature.id))

        for feature in other.features:
            new_record.features.append(SeqFeatureEM2(parent=new_record,
                                                     location=FeatureLocation(
                                                         feature.location.start + len(
                                                             self.seq) + offset,
                                                         feature.location.end + len(
                                                             self.seq) + offset),
                                                     strand=feature.strand, id=feature.id))

        return new_record

    def stitch(self, other, fpos_in_self, fpos_in_other, feature_length, **kwargs):
        """
        Stitches two records, that is, joins them according to an overlapping feature. The sequences
        may or may not overlap. If not, Ns or Xs are added to fill the gap. If they overlap, a
        warning is issued if sequences do not correspond exactly. The new record keeps track of the
        two original records as Features. By convention, the self record should contain the start
        position of the feature, the other contains the end position and the overlapping feature
        should be on the same strand in both records. It is the user's responsibility to provide
        the records in the right order and direction.
        :param other: the other SeqRecordEM2 object to stich
        :param fpos_in_self: feature position in self record (start position of feature)
        :param fpos_in_other: feature position in other record (end position of feature)
        :param feature_length: feature length
        :param kwargs: any additional parameters that may be passed to the created record
        :return: the stitched record as a new SeqRecordEM2 object
        """
        # TODO define overlapping feature as a new feature of the resulting stitched record
        # TODO define original records as features of the resulting stitched record
        offset = feature_length + fpos_in_self - fpos_in_other - len(self.seq) - 1
        print(fpos_in_self, fpos_in_other, feature_length, len(self.seq), offset)
        return self.join(other, offset)