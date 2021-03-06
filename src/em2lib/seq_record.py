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

    def add_feature(self, **kwargs):
        """
        Adds a feature to the current record according to arguments passed as \*\*kwargs.

        :param kwargs: keyword arguments to pass to SeqFeatureEM2 class
        """
        self.features.append(SeqFeatureEM2(parent=self, **kwargs))

    def overlap(self, start, end=None, strand=0):
        """
        Retrieves features that overlap a given position range.

        :param strand: strand specification of features to be returned. If strand is 0, then\
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
        strand or both. If nearest is True, then only the nearest ones are returned.

        :param nearest: if True, only the nearest features are returned. This only makes sense when\
        strand is 0
        :param strand: strand specification of features to be returned. If strand is 0, then\
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
        one strand or both. If nearest is True, then only the nearest ones are returned.

        :param nearest: if True, only the nearest features are returned. This only makes sense when\
        strand is 0
        :param strand: strand specification of features to be returned. If strand is 0, then\
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
        :param strand: strand specification of features to be returned. If strand is 0, then\
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

        :param keepself: if True and overlapping subsequences are different, then keep sequence\
        from self record, otherwise keep the sequence of other record.
        :param other: the other SeqRecordEM2 object
        :param offset: the offset of the two sequences. If the value is negative, then the two\
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

        new_record.id = self.id + '_' + other.id

        for feature in self.features:
            new_record.add_feature(location=feature.location, strand=feature.strand, id=feature.id)

        for feature in other.features:
            new_record.add_feature(
                location=FeatureLocation(feature.location.start + len(self.seq) + offset,
                                         feature.location.end + len(self.seq) + offset),
                strand=feature.strand, id=feature.id)
        return new_record

    def stitch(self, other, fpos_in_self, fpos_in_other, feature_length, orientation=1, **kwargs):
        """
        Stitches two records, that is, joins them according to an overlapping feature. The sequences
        may or may not overlap. If not, Ns or Xs are added to fill the gap. If they overlap, a
        warning is issued if sequences do not correspond exactly. The new record keeps track of the
        two original records as Features. By convention, the self record should contain the start
        position of the feature on the forward strand, the other contains the end position on the
        forward strand if orientation=1 or on the reverse strand if orientation=-1. If orientation
        is -1, then the other record is reversed/complemented before stitching and the position of
        the overlapping feature is modified accordingly. It is the user's responsibility to provide
        the records in the right order.

        :param other: the other SeqRecordEM2 object to stitch
        :param fpos_in_self: feature position in self record (start position of feature)
        :param fpos_in_other: feature position in other record (end position of feature), according
         to FeatureLocation conventions for end position requiring that length = end - start, it is
         not included in the feature.
        :param feature_length: feature length
        :param orientation: the orientation of the other record relative to the self, either 1 if
         it is in the same orientation, -1 if other needs to be reversed before stitching, 0 if
         stranded but unknown, None for proteins
        :param kwargs: any additional parameters that may be passed to the constructor of the
         stitching feature in the new record
        :return: the stitched record as a new SeqRecordEM2 object
        """
        orientation = None if other.seq.is_protein() else orientation
        if orientation == -1:
            other = other.reverse_complement()
            fpos_in_other = len(other.seq) - fpos_in_other - 1
        offset = feature_length + fpos_in_self - fpos_in_other - len(self.seq)
        stitched = self.join(other, offset)
        # Adding self sequence as a feature of new record
        stitched.add_feature(location=FeatureLocation(0, len(self.seq)), id=self.id)
        stitched.features[-1].strand = None if self.seq.is_protein() else 1

        # Adding other sequence as a feature of new record
        stitched.add_feature(location=FeatureLocation(len(self.seq) + offset,
                                                      len(self.seq) + len(other.seq) + offset),
                             id=other.id, strand=orientation)

        # Adding stitching feature as a feature of new record
        stitched.add_feature(location=FeatureLocation(fpos_in_self, fpos_in_self + feature_length),
                             **kwargs)
        if self.seq.is_protein():
            stitched.features[-1].strand = None
        return stitched

    def reverse_complement(self, id=False, name=False, description=False, features=True,
                           annotations=False, letter_annotations=True, dbxrefs=False):
        """
        Reverse-complement the record adjusting features and their positions accordingly. The record
        id is conserved but if name is not specified 'reversed' is appended. All other arguments are
        passed and handled by the parent method.
        Note that the main goal for this method is to replace SeqRecord and Seq objects by their
        SeqRecordEM2 and SeqEM2 equivalents when reverse/complementing.

        :param id: the id for the reversed record
        :param name: the name for the reversed record
        :param description: the description for the reversed record
        :param features: keep and adjust location of features if True
        :param annotations: keep annotations if True
        :param letter_annotations: keep letter_annotations if True
        :param dbxrefs: keep dbxrefs if True
        :return: a reversed copy of the record
        """
        id = self.id if id is False else id
        name = self.name + ' reversed' if name is False else name
        rev_record = super().reverse_complement(id=id, name=name, description=self.description,
                                                features=features, annotations=annotations,
                                                letter_annotations=letter_annotations,
                                                dbxrefs=dbxrefs)
        rev_record = SeqRecordEM2(SeqEM2.dna(str(rev_record.seq)), id=self.id, name=rev_record.name,
                                  description=self.description,
                                  features=rev_record.features, annotations=rev_record.annotations,
                                  letter_annotations=rev_record.letter_annotations,
                                  dbxrefs=rev_record.dbxrefs)
        return rev_record

    def orfs_to_features(self, start=['ATG'], stop=None, filter=None, add=False):
        """
        Determines all open reading frames in a sequence record. All the returned ORFs have a length
        that is a multiple of 3. Thus, for sequences without any stop codon, 3 ORFs are returned,
        one for each frame. Both strands are examined but it is possible to filter the ORFs by
        length, frame, etc. with the FeatureFilter defined by the filter argument.

        :param start: a list of accepted start codons
        :param stop: a list of accepted stop codons
        :param filter: a FeatureFilter object defining a filter to select ORFs according to some
         criteria
        :param add: if True, the selected ORFs are added to the record's features
        :return: a list of ORFs as SeqFeatureEM2 objects
        """
        features = []
        for orf in self.seq.get_orfs(start=start, stop=stop):
            features.append(
                SeqFeatureEM2(location=FeatureLocation(orf[0], orf[1]), type='ORF', strand=1))
        for orf in self.reverse_complement().seq.get_orfs(start=start, stop=stop):
            features.append(
                SeqFeatureEM2(location=FeatureLocation(orf[0], orf[1]), type='ORF', strand=-1))
        if filter is not None:
            features = filter.apply(features)
        if add:
            self.features.extend(features)
        return features
