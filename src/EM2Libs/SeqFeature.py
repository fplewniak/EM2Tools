"""
Extension of Bio.SeqFeature module from Biopython to add or improve functionalities
"""
#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
import warnings

from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition


class SeqFeatureEM2(SeqFeature):
    """
    SeqFeatureEM2 class providing extension to Bio.SeqFeature.SeqFeature class of BioPython package.
    """

    def __init__(self, parent=None, ref=None, **kwargs):
        super().__init__(**kwargs)
        self.parent = parent
        self.ref = parent.id if (parent is not None) and (ref is None) else ref

    def lies_within(self, start, end):
        """
        Determines whether feature lies entirely with the specified range. Fuzzy positions are turned into integers.

        :param start: start of range either int or ExactPosition
        :param end: end of range either int or ExactPosition

        :return: True if feature boundaries lie with the specified range.
        """
        if not all([isinstance(x, ExactPosition) for x in [start, end, self.location.start, self.location.end]]):
            warnings.warn('Warning: fuzzy positions are not handled as such but are converted into integer values.', )
        return start <= self.location.start <= self.location.end <= end

    def covers(self, start, end):
        """
        Determines whether feature covers the whole range specified by start and end

        :param start: start of range either int or ExactPosition
        :param end: end of range either int or ExactPosition, if None then end=start

        :return: True if feature covers the specified range
        """
        return self.overlaps(start) and self.overlaps(end)

    def overlaps(self, start, end=None):
        """
        Determines whether feayur overlaps a position range.

        :param start: start of range either int or ExactPosition
        :param end: end of range either int or ExactPosition

        :return: True if feature overlaps range
        """
        if end is None:
            end = start
        if not all([isinstance(x, ExactPosition) for x in [start, end, self.location.start, self.location.end]]):
            warnings.warn('Warning: fuzzy positions are not handled as such but are converted into integer values.', )
        left = max(self.location.start, start)
        right = min(self.location.end, end)
        return left <= right

    def intersect(self, other, **kwargs):
        """
        Creates a new feature which is the intersection of feature and another one

        :param other: the other feature
        """
        if id(self.parent) != id(other.parent):
            raise ValueError('Undetermined intersection of features because of different parent sequence records')
        if self.overlaps(other.location.start, other.location.end):
            start = max(self.location.start, other.location.start)
            end = min(self.location.end, other.location.end)
            return SeqFeatureEM2(parent=self.parent, location=FeatureLocation(start, end), **kwargs)
        return None
