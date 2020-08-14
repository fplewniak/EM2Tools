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

    def __init__(self, parent=None, **kwargs):
        super().__init__(**kwargs)
        self.parent = parent

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
        Determines whether feature lies entirely with the specified range. Fuzzy positions are turned into integers.

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


class FeatureFilter:
    """
    A class for the definition and application of a filter to a list of features. May be used to
    retrieve features from a record or any list of features according to length, location, type,
    strand and frame. It is possible to combine several criteria into one single filter.
    """
    def __init__(self):
        self._minlength = None
        self._maxlength = None
        self._strand = None
        self._frame = None
        self._covers = None
        self._overlaps = None
        self._lies_within = None
        self._keep = True
        self._type = None

    def keep(self, keep=True):
        """
        Set the keep attribute of the filter. If True, the features corresponding to the specified
        criteria will be kept, otherwise, they will be discarded.

        :param keep: boolean, True if features consistent with the criteria should be kept.
        :return: the current filter
        """
        self._keep = keep
        return self

    def type(self, feat_type=None):
        """
        Set the type of feature.

        :param feat_type: the type of feature to select
        :return: the current filter
        """
        self._type = feat_type
        return self

    def length(self, minlength=None, maxlength=None):
        """
        Set the length range of the feature.

        :param minlength: the minimum length of the feature or 0 if None
        :param maxlength: the maximum length of the feature or no limit if set to None
        :return: the current filter
        """
        self._minlength = minlength
        self._maxlength = maxlength
        return self

    def strand(self, strand=0):
        """
        Set the strand of the feature.

        :param strand: the strand of the feature
        :return: the current filter
        """
        self._strand = strand
        return self

    def frame(self, frame=0, strand=1):
        """
        Set the frame of the filter. By default, this is on the forward strand unless the strand is set to -1.

        :param frame: the frame 0, 1 or 2
        :param strand: the strand 1 or - 1
        :return: the current filter
        """
        self._frame = frame
        self._strand = strand
        return self

    def covers(self, start=None, end=None):
        """
        Set a region that must be covered by the feature.

        :param start: start position of the region
        :param end: end position of the region
        :return: the current filter
        """
        self._covers = (start, end)
        return self

    def overlaps(self, start=None, end=None):
        """
        Set a region that must overlap the feature.

        :param start: start position of the region
        :param end: end position of the region
        :return: the current filter
        """
        self._overlaps = (start, end)
        return self

    def lies_within(self, start=None, end=None):
        """
        Set a region within which the feature must lie.

        :param start: start position of the region
        :param end: end position of the region
        :return: the current filter
        """
        self._lies_within = (start, end)
        return self

    def type_applies(self, feature):
        """
        Test if type criterion applies to the feature and return a boolean stating whether the
        feature should be kept or rejected. If the criterion has not been set, then True is returned.

        :param feature: the feature to test
        :return: True the feature must be kept
        """
        if self._type is None:
            return True
        if self._keep is True:
            return feature.type == self._type
        return not feature.type == self._type

    def length_applies(self, feature):
        """
        Test if length criterion applies to the feature and return a boolean stating whether the
        feature should be kept or rejected. If the criterion has not been set, then True is returned.

        :param feature: the feature to test
        :return: True the feature must be kept
        """
        if self._keep is True:
            return feature.length_in_range(self._minlength, self._maxlength)
        return (self._minlength, self._maxlength) == (None, None) or feature.length_in_range(
            self._minlength,
            self._maxlength) is False

    def covers_applies(self, feature):
        """
        Test if covers criterion applies to the feature and return a boolean stating whether the
        feature should be kept or rejected. If the criterion has not been set, then True is returned.

        :param feature: the feature to test
        :return: True the feature must be kept
        """
        if self._covers is None:
            return True
        if self._keep is True:
            return feature.covers(self._covers[0], self._covers[1])
        return not feature.covers(self._covers[0], self._covers[1])

    def overlaps_applies(self, feature):
        """
        Test if overlaps criterion applies to the feature and return a boolean stating whether the
        feature should be kept or rejected. If the criterion has not been set, then True is returned.

        :param feature: the feature to test
        :return: True the feature must be kept
        """
        if self._overlaps is None:
            return True
        if self._keep is True:
            return feature.overlaps(self._overlaps[0], self._overlaps[1])
        return not feature.overlaps(self._overlaps[0], self._overlaps[1])

    def lies_within_applies(self, feature):
        """
        Test if lies_within criterion applies to the feature and return a boolean stating
        the feature should be kept or rejected. If the criterion has not been set, then True is returned.

        :param feature: the feature to test
        :return: True the feature must be kept
        """
        if self._lies_within is None:
            return True
        if self._keep is True:
            return feature.lies_within(self._lies_within[0], self._lies_within[1])
        return not feature.lies_within(self._lies_within[0], self._lies_within[1])

    def location_applies(self, feature):
        """
        Test if location criterion applies to the feature and return a boolean stating whether the
        feature should be kept or rejected. If the criterion has not been set, then True is returned.

        :param feature: the feature to test
        :return: True the feature must be kept
        """
        return all([self.covers_applies(feature),
                    self.overlaps_applies(feature),
                    self.lies_within_applies(feature)])

    def strand_applies(self, feature):
        """
        Test if strand criterion applies to the feature and return a boolean stating whether the
        feature should be kept or rejected. If the criterion has not been set, then True is
        returned.

        :param feature: the feature to test
        :return: True the feature must be kept
        """
        if self._strand is None:
            return True
        if self._keep is True:
            return feature.strand == self._strand
        return feature.strand != self._strand

    def frame_applies(self, feature):
        """
        Test if frame criterion applies to the feature and return a boolean stating whether the
        feature should be kept or rejected. If the criterion has not been set, then True is
        returned.

        :param feature: the feature to test
        :return: True the feature must be kept
        """
        if self._frame is None:
            return True
        pos = feature.location.start if self._strand == 1 else feature.location.end
        if self._keep is True:
            return (self._frame == pos % 3) and (feature.strand == self._strand)
        return not (self._frame == pos % 3 and feature.strand == self._strand)

    def apply(self, features):
        """
        Test if all defined criteria apply to the features in the list and return the list of
        features corresponding to the specified criteria.

        :param features: the list of features to filter
        :return: the filtered list of features
        """
        filtered = []
        for feat in features:
            if all([self.length_applies(feat),
                    self.strand_applies(feat),
                    self.location_applies(feat),
                    self.type_applies(feat),
                    self.frame_applies(feat)
                    ]):
                filtered.append(feat)
        return filtered
