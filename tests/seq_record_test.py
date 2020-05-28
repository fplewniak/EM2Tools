#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
import pytest


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


def test_join_sequences(dna_rec1, dna_rec2, dna_rec3):
    assert str(dna_rec1.join(dna_rec2, offset=3).seq) == str(dna_rec1.seq) + 'NNN' + str(dna_rec2.seq)
    assert str(dna_rec2.join(dna_rec1, offset=-73).seq) == str(dna_rec1.seq) + 'NNN' + str(dna_rec2.seq)
    assert str(
        dna_rec1.join(dna_rec2, offset=-10).seq) == 'ATGAGTCGGTAACGATGCATGCATGCAGCTGACGCATGAGTCGGTAACGATGCATGCATG'
    assert str(
        dna_rec2.join(dna_rec1, offset=-60).seq) == 'ATGAGTCGGTAACGATGCATGCATGCAGCTGACGCATGAGTCGGTAACGATGCATGCATG'
    assert str(
        dna_rec1.join(dna_rec3, offset=-10).seq) == 'ATGAGTCGGTAACGATGCATGCATGCAGCTGACGCATGAGTCGGTAACGATGCATGCATG'
    assert str(dna_rec1.join(dna_rec3, offset=-10,
                             keepself=False).seq) == 'ATGAGTCGGTAACGATGCATGCATGCACCTGACGCATGAGTCGGTAACGATGCATGCATG'
    assert str(
        dna_rec3.join(dna_rec1, offset=-60).seq) == 'ATGAGTCGGTAACGATGCATGCATGCACCTGACGCATGAGTCGGTAACGATGCATGCATG'
    assert str(dna_rec3.join(dna_rec1, offset=-60,
                             keepself=False).seq) == 'ATGAGTCGGTAACGATGCATGCATGCAGCTGACGCATGAGTCGGTAACGATGCATGCATG'


def test_join_differences(dna_rec1, dna_rec2, dna_rec3, protein_rec):
    with pytest.warns(UserWarning, match=r'.*Overlapping subsequences are different.*'):
        assert dna_rec1.join(dna_rec3, offset=-10)
    with pytest.raises(ValueError, match=r'Sequences are not of the same type.*'):
        assert dna_rec1.join(protein_rec)


def test_join_features(dna_rec1, dna_rec2, dna_rec3):
    newrecord = dna_rec1.join(dna_rec1, offset=3)
    assert [(f.id, int(f.location.start), int(f.location.end), f.strand, f.ref) for f in newrecord.features] == [
        ('A1', 0, 5, 1, '<unknown id>'),
        ('B1', 28, 33, 1, '<unknown id>'),
        ('A1', 38, 43, 1, '<unknown id>'),
        ('B1', 66, 71, 1, '<unknown id>'),
    ]

    newrecord = dna_rec1.join(dna_rec1, offset=-73)
    assert [(f.id, int(f.location.start), int(f.location.end), f.strand, f.ref) for f in newrecord.features] == [
        ('A1', 0, 5, 1, '<unknown id>'),
        ('B1', 28, 33, 1, '<unknown id>'),
        ('A1', 38, 43, 1, '<unknown id>'),
        ('B1', 66, 71, 1, '<unknown id>'),
    ]

    newrecord = dna_rec1.join(dna_rec2, offset=-10)
    assert [(f.id, int(f.location.start), int(f.location.end), f.strand, f.ref) for f in newrecord.features] == [
        ('A1', 0, 5, 1, '<unknown id>'),
        ('B1', 28, 33, 1, '<unknown id>'),
        ('A2', 25, 30, 1, '<unknown id>'),
        ('B2', 53, 58, 1, '<unknown id>'),
    ]

    newrecord = dna_rec2.join(dna_rec1, offset=-60)
    assert [(f.id, int(f.location.start), int(f.location.end), f.strand, f.ref) for f in newrecord.features] == [
        ('A1', 0, 5, 1, '<unknown id>'),
        ('B1', 28, 33, 1, '<unknown id>'),
        ('A2', 25, 30, 1, '<unknown id>'),
        ('B2', 53, 58, 1, '<unknown id>'),
    ]


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
