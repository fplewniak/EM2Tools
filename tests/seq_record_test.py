#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
import warnings

import pytest
from em2lib.seq_feature import FeatureFilter as FF, SeqFeatureEM2
from Bio.SeqFeature import FeatureLocation


def test_overlap(dna_rec):
    assert sorted([f.id for f in dna_rec.overlap(16, 25)]) == ['C', 'D', 'G', 'H']
    assert sorted([f.id for f in dna_rec.overlap(16)]) == ['C', 'H']
    assert sorted([f.id for f in dna_rec.overlap(16, 25, strand=-1)]) == ['G', 'H']
    assert sorted([f.id for f in dna_rec.overlap(12, 14, strand=-1)]) == []


def test_features_after(dna_rec):
    assert sorted([f.id for f in dna_rec.feature_after(16)]) == ['D', 'F']
    assert sorted([f.id for f in dna_rec.feature_after(16, nearest=True)]) == ['D']
    assert sorted([f.id for f in dna_rec.feature_after(16, strand=1)]) == ['D']
    assert sorted([f.id for f in dna_rec.feature_after(16, strand=-1)]) == ['F']
    assert sorted([f.id for f in dna_rec.feature_after(14)]) == ['C', 'F']
    assert sorted([f.id for f in dna_rec.feature_after(14, nearest=True)]) == ['C']
    assert sorted([f.id for f in dna_rec.feature_after(14, strand=1)]) == ['C']
    assert sorted([f.id for f in dna_rec.feature_after(14, strand=-1)]) == ['F']
    assert sorted([f.id for f in dna_rec.feature_after(15)]) == ['F', 'H']
    assert sorted([f.id for f in dna_rec.feature_after(15, nearest=True)]) == ['H']
    assert sorted([f.id for f in dna_rec.feature_after(15, strand=1)]) == ['H']
    assert sorted([f.id for f in dna_rec.feature_after(15, strand=-1)]) == ['F']
    assert sorted([f.id for f in dna_rec.feature_after(4, strand=-1)]) == []


def test_features_before(dna_rec):
    assert sorted([f.id for f in dna_rec.feature_before(16)]) == ['B', 'G']
    assert sorted([f.id for f in dna_rec.feature_before(16, nearest=True)]) == ['G']
    assert sorted([f.id for f in dna_rec.feature_before(16, strand=1)]) == ['B']
    assert sorted([f.id for f in dna_rec.feature_before(16, strand=-1)]) == ['G']
    assert sorted([f.id for f in dna_rec.feature_before(20)]) == ['H']
    assert sorted([f.id for f in dna_rec.feature_before(20, strand=1)]) == ['H']
    assert sorted([f.id for f in dna_rec.feature_before(20, strand=-1)]) == []


def test_surrounding(dna_rec):
    assert sorted([f.id for f in dna_rec.surrounding_features(19)]) == ['B', 'D', 'F']
    assert sorted([f.id for f in dna_rec.surrounding_features(19, nearest=True)]) == ['D']
    assert sorted([f.id for f in dna_rec.surrounding_features(16)]) == ['B', 'D', 'F', 'G']
    assert sorted([f.id for f in dna_rec.surrounding_features(16, strand=1)]) == ['B', 'D']
    assert sorted([f.id for f in dna_rec.surrounding_features(16, strand=-1)]) == ['F', 'G']
    assert sorted([f.id for f in dna_rec.surrounding_features(15)]) == ['B', 'F', 'H']
    assert sorted([f.id for f in dna_rec.surrounding_features(4)]) == ['A', 'B', 'F']
    assert sorted([f.id for f in dna_rec.surrounding_features(26)]) == ['C', 'G']
    assert sorted([f.id for f in dna_rec.surrounding_features(26, nearest=True)]) == ['G']
    assert sorted([f.id for f in dna_rec.surrounding_features(16, nearest=True)]) == ['G']
    assert sorted([f.id for f in dna_rec.surrounding_features(16, nearest=True, strand=1)]) == ['D']
    assert sorted([f.id for f in dna_rec.surrounding_features(16, nearest=True, strand=-1)]) == ['G']
    assert sorted([f.id for f in dna_rec.surrounding_features(11, nearest=True)]) == ['B', 'F']
    assert sorted([f.id for f in dna_rec.surrounding_features(11, nearest=True, strand=1)]) == ['B']
    assert sorted([f.id for f in dna_rec.surrounding_features(11, nearest=True, strand=-1)]) == ['F']


def test_add_feature(dna_rec1_NNN_rec2, dna_stitch_rec1_NNN_rec2):
    dna_rec1_NNN_rec2.add_feature(location=FeatureLocation(0, 35), strand=1, id='Rec1')
    dna_rec1_NNN_rec2.add_feature(location=FeatureLocation(38, 73), strand=1, id='Rec2')
    dna_rec1_NNN_rec2.add_feature(location=FeatureLocation(30, 48), strand=1, id='stitcher')
    assert [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in dna_rec1_NNN_rec2.features] == \
           [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in dna_stitch_rec1_NNN_rec2.features]


def test_join_sequences(dna_rec1, dna_rec2, dna_rec3, dna_rec1_NNN_rec2, dna_rec1_overlap_rec2,
                        dna_rec1_overlap_rec3, dna_rec1_overlap_rec3_keep_false):
    assert str(dna_rec1.join(dna_rec2, offset=3).seq) == str(dna_rec1_NNN_rec2.seq)
    assert str(dna_rec2.join(dna_rec1, offset=-73).seq) == str(dna_rec1_NNN_rec2.seq)
    assert str(dna_rec1.join(dna_rec2, offset=-10).seq) == str(dna_rec1_overlap_rec2.seq)
    assert str(dna_rec2.join(dna_rec1, offset=-60).seq) == str(dna_rec1_overlap_rec2.seq)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert str(dna_rec1.join(dna_rec3, offset=-10).seq) == str(dna_rec1_overlap_rec3.seq)
        assert str(dna_rec1.join(dna_rec3, offset=-10, keepself=False).seq) == str(dna_rec1_overlap_rec3_keep_false.seq)
        assert str(dna_rec3.join(dna_rec1, offset=-60).seq) == str(dna_rec1_overlap_rec3_keep_false.seq)
        assert str(dna_rec3.join(dna_rec1, offset=-60, keepself=False).seq) == str(dna_rec1_overlap_rec3.seq)


def test_join_differences(dna_rec1, dna_rec2, dna_rec3, protein_rec):
    with pytest.warns(UserWarning, match=r'.*Overlapping subsequences are different.*'):
        assert dna_rec1.join(dna_rec3, offset=-10)
    with pytest.raises(ValueError, match=r'Sequences are not of the same type.*'):
        assert dna_rec1.join(protein_rec)


def test_join_features(dna_rec1, dna_rec2, dna_rec3, dna_rec1_NNN_rec2, dna_rec1_overlap_rec2):
    newrecord = dna_rec1.join(dna_rec2, offset=3)
    assert [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in newrecord.features] == \
           [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in dna_rec1_NNN_rec2.features]

    newrecord = dna_rec2.join(dna_rec1, offset=-73)
    assert [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in newrecord.features] == \
           [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in dna_rec1_NNN_rec2.features]

    newrecord = dna_rec1.join(dna_rec2, offset=-10)
    assert [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in newrecord.features] == \
           [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in dna_rec1_overlap_rec2.features]

    newrecord = dna_rec2.join(dna_rec1, offset=-60)
    assert [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in newrecord.features] == \
           [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in dna_rec1_overlap_rec2.features]


def test_stitch_sequences(dna_rec1, dna_rec2, dna_rec1_overlap_rec2, dna_rec1_NNN_rec2):
    assert str(dna_rec1.stitch(dna_rec2, 22, 20, 23).seq) == str(dna_rec1_overlap_rec2.seq)
    assert str(dna_rec1.stitch(dna_rec2, 30, 10, 18).seq) == str(dna_rec1_NNN_rec2.seq)


def test_stitch_features(dna_rec1, dna_rec2, dna_rec2_rev, dna_stitch_rec1_overlap_rec2,
                         dna_stitch_rec1_NNN_rec2, dna_stitch_rec1_NNN_rec2rev):
    newrecord = dna_rec1.stitch(dna_rec2, 22, 20, 23, orientation=1, strand=1, id='stitcher')
    assert [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in newrecord.features] == \
           [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in dna_stitch_rec1_overlap_rec2.features]

    newrecord = dna_rec1.stitch(dna_rec2, 30, 10, 18, orientation=1, strand=1, id='stitcher')
    assert [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in newrecord.features] == \
           [(f.id, int(f.location.start), int(f.location.end), f.strand) for f in dna_stitch_rec1_NNN_rec2.features]

    newrecord = dna_rec1.stitch(dna_rec2_rev, 30, 24, 18, orientation=-1, strand=1, id='stitcher')
    assert [(f.id, int(f.location.start), int(f.location.end), f.strand, f.ref) for f in newrecord.features] == \
           [(f.id, int(f.location.start), int(f.location.end), f.strand, f.ref) for f in dna_stitch_rec1_NNN_rec2rev.features]


def test_comparisons(dna_rec1, dna_rec2, dna_rec3):
    assert dna_rec1.__lt__(dna_rec2)
    assert not dna_rec3.__lt__(dna_rec2)
    assert dna_rec1.__le__(dna_rec2)
    assert dna_rec1.__le__(dna_rec1)
    assert not dna_rec2.__le__(dna_rec1)
    assert dna_rec1.__eq__(dna_rec1)
    assert not dna_rec3.__eq__(dna_rec2)
    assert dna_rec1.__ne__(dna_rec2)
    assert not dna_rec3.__ne__(dna_rec3)
    assert dna_rec2.__gt__(dna_rec1)
    assert not dna_rec1.__gt__(dna_rec2)
    assert dna_rec3.__ge__(dna_rec2)
    assert dna_rec1.__ge__(dna_rec1)
    assert not dna_rec1.__ge__(dna_rec2)
    assert dna_rec3 >= dna_rec3 > dna_rec1 == dna_rec1 < dna_rec2 <= dna_rec2
    assert sorted([dna_rec3, dna_rec1, dna_rec2]) == [dna_rec1, dna_rec2, dna_rec3]


def test_orfs_to_features(dna_rec1_overlap_rec2):
    assert {(int(f.location.start), int(f.location.end), f.strand, f.type)
            for f in dna_rec1_overlap_rec2.orfs_to_features(filter=FF().length(30))} == \
           {(22, 58, 1, 'ORF'), (1, 58, -1, 'ORF'), (5, 59, -1, 'ORF'), (23, 59, -1, 'ORF')}
    assert {(int(f.location.start), int(f.location.end), f.strand, f.type)
            for f in dna_rec1_overlap_rec2.orfs_to_features(start=None,
                                                            filter=FF().length(30))} == \
           {(4, 58, 1, 'ORF'), (1, 58, -1, 'ORF'), (2, 59, -1, 'ORF'), (0, 60, -1, 'ORF')}
    dna_rec1_overlap_rec2.orfs_to_features(filter=FF().length(30), add=True)
    assert {(int(f.location.start), int(f.location.end), f.strand, f.type)
            for f in dna_rec1_overlap_rec2.features} == \
           {(22, 58, 1, 'ORF'), (1, 58, -1, 'ORF'), (5, 59, -1, 'ORF'), (23, 59, -1, 'ORF'),
            (0, 5, 1, ''), (28, 33, 1, ''), (25, 30, 1, ''), (53, 58, 1, '')}