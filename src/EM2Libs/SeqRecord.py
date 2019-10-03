"""
Extension module to the Biopython Bio.SeqRecord module
"""

#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
from Bio.SeqRecord import SeqRecord


class SeqRecordEM2(SeqRecord):
    """
    Extension to Biopython SeqRecord class
    """

    def __init__(self, seq, **kwargs):
        super().__init__(seq, **kwargs)

    def overlap(self, start, end=None, strand=0):
        """
        Retrieves features that overlap a given position range.
        :param strand: strand specification of features to be returned. If strand is 0, then features on
        both strands are returned
        :param start: start of range
        :param end: end of range, if None, then end=start
        """
        if end is None:
            end = start
        if strand == 0:
            return [f for f in self.features if f.overlaps(start, end)]
        return [f for f in self.features if f.overlaps(start, end) and f.strand == strand]

    def feature_after(self, position, strand=1, nearest=False):
        """
        Retrieves the features immediately after (but not overlaping) the specified position, on one strand or both.
        If nearest is True, then only the nearest ones are returned.
        :param nearest:
        :param strand: strand specification of features to be returned. If strand is 0, then features on
         both strands are returned
        :param position: the position
        """
        raise NotImplementedError

    def feature_before(self, position, strand=1, nearest=False):
        """
        Retrieves the features immediately before (but not overlaping) the specified position, on one strand or both.
        If nearest is True, then only the nearest ones are returned.
        :param nearest: if True, only the nearest features are returned. This only makes sense when strand is O
        :param strand: strand specification of features to be returned. If strand is 0, then features on
         both strands are returned
        :param position: the position
        """
        raise NotImplementedError

    def surrounding_features(self, position, strand=0, nearest=False):
        """
        Retrieves all the features around a given position but not overlaping it. If nearest is True, then only the
        nearest features are returned.
        :param position: the position
        :param nearest: if True, only the nearest features are returned. This only makes sense when strand is O
        :param strand: strand specification of features to be returned. If strand is 0, then features on
         both strands are returned
        """
        raise NotImplementedError

    def stitch(self, other=None, offset=0):
        """
        Joins two SeqRecordEM2 objects into a new one representing the resulting merged sequence
        :param other: the other SeqRecordEM2 object
        :param offset: the offset of the two sequences. If the value is negative, then the two sequences overlap.
        """
        raise NotImplementedError
