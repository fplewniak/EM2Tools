"""
Extension of Bio.SeqFeature module from Biopython to add or improve functionalities
"""
#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#

from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition, BeforePosition, AfterPosition


class SeqFeatureEM2(SeqFeature):
    """
    SeqFeatureEM2 class providing extension to Bio.SeqFeature.SeqFeature class of BioPython package.
    """

    def __init__(self, parent=None, **kwargs):
        super().__init__(**kwargs)
        self.parent = parent

    def lies_within(self, start, end):
        """
        Determines whether feature lies entirely with the specified range
        :param start: start of range either int, ExactPosition, BeforePosition or AfterPosition
        :param end: end of range either int, ExactPosition, BeforePosition or AfterPosition
        :return: the probability that feature lies within the specified range
        """
        if isinstance(start, int):
            start = ExactPosition(start)
        if isinstance(end, int):
            end = ExactPosition(end)
        if all([isinstance(p, ExactPosition) for p in [start, end, self.location.start, self.location.end]]):
            return start <= self.location.start <= self.location.end <= end
        # return start <= self.location.nofuzzy_start <= self.location.nofuzzy_end <= end
        raise NotImplementedError('Not implemented')

    def overlaps(self, position):
        """
        Determines whether feature overlaps a given position
        :param position: the position either int, ExactPosition, BeforePosition or AfterPosition
        :return: the probability that feature overlaps the specified position
        """
        print(position, self.location.start, self.location.end)
        if isinstance(position, int):
            position = ExactPosition(position)
        if self.location.nofuzzy_start <= position <= self.location.nofuzzy_end:
            return 1.0
        if all([isinstance(p, ExactPosition) for p in [position, self.location.start, self.location.end]]) is False:
            if isinstance(self.location.start, BeforePosition) and position < self.location.nofuzzy_start:
                return 1.0 / self.location.nofuzzy_start
            if isinstance(self.location.end, AfterPosition) and position > self.location.nofuzzy_end:
                return 1.0 / self.location.nofuzzy_end
        return 0.0

    def covers(self, start, end):
        """
        Determines whether feature covers the whole range specified by start and end
        :param start: start of range either int, ExactPosition, BeforePosition or AfterPosition
        :param end: end of range either int, ExactPosition, BeforePosition or AfterPosition
        :return: the probability that feature covers the specified range
        """
        return self.overlaps(start) * self.overlaps(end)

    def intersect(self, other, *args, **kwargs):
        """
        Creates a new feature which is the intersection of feature and another one
        :param other: the other feature
        :param ftype: the type for the newly created feature
        :param fid: the id for the newly created feature
        """
        if id(self.parent) != id(other.parent):
            raise ValueError('Undetermined intersection of features because of different parent sequence records')
        left, right = [self, other] if self.location.nofuzzy_start <= other.location.nofuzzy_start else [other, self]
        print(left.location, right.location, left.overlaps(right.location.start), right.overlaps(left.location.end))
        if 0.0 < left.overlaps(right.location.start) < 1.0 or 0.0 < right.overlaps(left.location.end) < 1.0:
            raise ValueError('Undetermined intersection of features due to fuzzy boundaries')

        return SeqFeatureEM2(parent=self.parent,
                             location=FeatureLocation(right.location.start, min(left.location.end, right.location.end)),
                             *args, **kwargs)
