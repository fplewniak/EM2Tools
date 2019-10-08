"""
Extension module to the Biopython Bio.SeqRecord module
"""

#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from EM2Libs.Seq import SeqEM2
from EM2Libs.SeqFeature import SeqFeatureEM2


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
        :return: a list of overlaping features
        :param strand: strand specification of features to be returned. If strand is 0, then features on
        both strands are returned. If feature.strand is 0, then all strands will match.
        :param start: start of range
        :param end: end of range, if None, then end=start
        """
        if end is None:
            end = start
        return [f for f in self.features if f.overlaps(start, end) and f.strand * strand >= 0]

    def feature_after(self, position, strand=0, nearest=False):
        """
        Retrieves the features immediately after (but not overlaping) the specified position, on one strand or both.
        If nearest is True, then only the nearest ones are returned.
        :param nearest: if True, only the nearest features are returned. This only makes sense when strand is O
        :param strand: strand specification of features to be returned. If strand is 0, then features on
         both strands are returned
        :param position: the position
        """
        after_fwd = {f: f.location.start for f in self.features if f.location.start > position and f.strand in [0, 1]}
        after_rev = {f: f.location.end for f in self.features if f.location.end < position and f.strand in [0, -1]}
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
        Retrieves the features immediately before (but not overlaping) the specified position, on one strand or both.
        If nearest is True, then only the nearest ones are returned.
        :param nearest: if True, only the nearest features are returned. This only makes sense when strand is O
        :param strand: strand specification of features to be returned. If strand is 0, then features on
         both strands are returned
        :param position: the position
        """
        before_fwd = {f: f.location.end for f in self.features if f.location.end < position and f.strand in [0, 1]}
        before_rev = {f: f.location.start for f in self.features if f.location.start > position and f.strand in [0, -1]}
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
        Retrieves all the features around a given position but not overlaping it. If nearest is True, then only the
        nearest features are returned.
        :param position: the position
        :param nearest: if True, only the nearest features are returned.
        :param strand: strand specification of features to be returned. If strand is 0, then features on
         both strands are returned
        """
        if nearest is False:
            return list(set(self.feature_after(position, strand) + self.feature_before(position, strand)))
        feat_list = {f: min(abs(position - f.location.start), abs(position - f.location.end)) for f in list(
            set(self.feature_after(position, strand, True) + self.feature_before(position, strand, True)))}
        return [f for f, v in feat_list.items() if v == min(feat_list.values())]

    def stitch(self, other=None, offset=0):
        """
        Joins two SeqRecordEM2 objects into a new one representing the resulting merged sequence
        :param other: the other SeqRecordEM2 object
        :param offset: the offset of the two sequences. If the value is negative, then the two sequences overlap.
        """

        if offset >= 0:
            if self.seq.is_protein():
                newseq = str(self.seq) + 'X' * offset + str(other.seq)
            else:
                newseq = str(self.seq) + 'N' * offset + str(other.seq)
        else:
            if str(self.seq)[offset:] != str(other.seq)[0:-offset]:
                raise ValueError('Warning!!! Overlapping subsequences are different.')
            newseq = str(self.seq) + str(other.seq)[-offset:]

        if self.seq.is_protein() and other.seq.is_protein():
            newrecord = SeqRecordEM2(SeqEM2.protein(newseq))
        elif not (self.seq.is_protein() or other.seq.is_protein()):
            newrecord = SeqRecordEM2(SeqEM2.dna(newseq))
        else:
            raise ValueError('Sequences are not of the same type. It is impossible to stitch them.')

        for f in self.features:
            newrecord.features.append(SeqFeatureEM2(parent=newrecord, location=f.location, strand=f.strand, id=f.id))

        for f in other.features:
            newrecord.features.append(SeqFeatureEM2(parent=newrecord,
                                                    location=FeatureLocation(f.location.start + len(self.seq) + offset,
                                                                             f.location.end + len(self.seq) + offset),
                                                    strand=f.strand, id=f.id))

        return newrecord
