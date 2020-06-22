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

    def length_in_range(self, minlength=None, maxlength=None):
        """
        Checks whether the feature length is with the specified range.

        :param minlength: lower length bound
        :param maxlength: upper length bound
        :return: True if feature length is within specified range, False otherwise
        """
        if minlength is None and maxlength is None:
            return True
        if minlength is None:
            return self.__len__() <= maxlength
        if maxlength is None:
            return minlength <= self.__len__()
        return minlength <= self.__len__() <= maxlength

    def lies_within(self, start, end):
        """
        Determines whether feature lies entirely with the specified range. Fuzzy positions are
         turned into integers.

        :param start: start of range either int or ExactPosition
        :param end: end of range either int or ExactPosition

        :return: True if feature boundaries lie with the specified range.
        """
        if not all([isinstance(x, ExactPosition) for x in
                    [start, end, self.location.start, self.location.end]]):
            warnings.warn(
                'Warning: fuzzy positions are not handled as such but are '
                'converted into integer values.', )
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
        Determines whether feature overlaps a position range.

        :param start: start of range either int or ExactPosition
        :param end: end of range either int or ExactPosition

        :return: True if feature overlaps range
        """
        if end is None:
            end = start
        if not all([isinstance(x, ExactPosition) for x in
                    [start, end, self.location.start, self.location.end]]):
            warnings.warn(
                'Warning: fuzzy positions are not handled as such but are'
                ' converted into integer values.', )
        left = max(self.location.start, start)
        right = min(self.location.end, end)
        return left <= right

    def intersect(self, other, **kwargs):
        """
        Creates a new feature which is the intersection of feature and another one

        :param other: the other feature
        """
        if id(self.parent) != id(other.parent):
            raise ValueError(
                'Undetermined intersection of features because of different'
                ' parent sequence records')
        if self.overlaps(other.location.start, other.location.end):
            start = max(self.location.start, other.location.start)
            end = min(self.location.end, other.location.end)
            return SeqFeatureEM2(parent=self.parent, location=FeatureLocation(start, end), **kwargs)
        return None

    def move(self, offset):
        """
        Moves a feature by a certain offset

        :param offset: offset by which the feature must be moved
        """
        self.location = FeatureLocation(self.location.start + offset, self.location.end + offset)
        return self


class FeatureFilter():
    def __init__(self):
        self._minlength = None
        self._maxlength = None
        self._strand = None
        self._covers = None
        self._overlaps = None
        self._lies_within = None
        self._keep = True

    def keep(self, keep=True):
        self._keep = keep
        return self

    def length(self, minlength=None, maxlength=None):
        self._minlength = minlength
        self._maxlength = maxlength
        return self

    def strand(self, strand=0):
        self._strand = strand
        return self

    def covers(self, start=None, end=None):
        self._covers = (start, end)
        return self

    def overlaps(self, start=None, end=None):
        self._overlaps = (start, end)
        return self

    def lies_within(self, start=None, end=None):
        self._lies_within = (start, end)
        return self

    def length_applies(self, feature):
        if self._keep is True:
            return feature.length_in_range(self._minlength, self._maxlength)
        return (self._minlength, self._maxlength) == (None, None) or feature.length_in_range(
            self._minlength,
            self._maxlength) is False

    def covers_applies(self, feature):
        pass

    def overlaps_applies(self, feature):
        pass

    def lies_within_applies(self, feature):
        pass

    def location_applies(self, feature):
        pass

    def strand_applies(self, feature):
        if self._strand is None:
            return True
        if self._keep is True:
            return feature.strand == self._strand
        return feature.strand != self._strand

    def apply(self, features):
        filtered = []
        for feat in features:
            if all([self.length_applies(feat),
                    self.strand_applies(feat)
                    ]):
                filtered.append(feat)
        return filtered
