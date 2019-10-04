#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#

import unittest

import pytest
from Bio.SeqFeature import FeatureLocation

from EM2Libs.Seq import SeqEM2
from EM2Libs.SeqFeature import SeqFeatureEM2
from EM2Libs.SeqRecord import SeqRecordEM2


class SeqRecordFeatureTests(unittest.TestCase):
    rec: SeqRecordEM2 = SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGC'), id='Rec', name='DummyDNA')

    rec.features = [
        SeqFeatureEM2(parent=rec, location=FeatureLocation(0, 2), strand=1, id='A'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(6, 10), strand=1, id='B'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(15, 20), strand=1, id='C'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(20, 30), strand=1, id='D'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(6, 10), strand=-1, id='F'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(18, 25), strand=-1, id='G'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(16, 19), strand=0, id='H')
    ]

    @classmethod
    def test_overlap(cls):
        assert sorted([f.id for f in cls.rec.overlap(16, 25)]) == ['C', 'D', 'G', 'H']
        assert sorted([f.id for f in cls.rec.overlap(16)]) == ['C', 'H']
        assert sorted([f.id for f in cls.rec.overlap(16, 25, strand=-1)]) == ['G', 'H']
        assert sorted([f.id for f in cls.rec.overlap(12, 14, strand=-1)]) == []

    @classmethod
    def test_features_after(cls):
        assert sorted([f.id for f in cls.rec.feature_after(16)]) == ['D', 'F']
        assert sorted([f.id for f in cls.rec.feature_after(16, nearest=True)]) == ['D']
        assert sorted([f.id for f in cls.rec.feature_after(16, strand=1)]) == ['D']
        assert sorted([f.id for f in cls.rec.feature_after(16, strand=-1)]) == ['F']
        assert sorted([f.id for f in cls.rec.feature_after(14)]) == ['C', 'F']
        assert sorted([f.id for f in cls.rec.feature_after(14, nearest=True)]) == ['C']
        assert sorted([f.id for f in cls.rec.feature_after(14, strand=1)]) == ['C']
        assert sorted([f.id for f in cls.rec.feature_after(14, strand=-1)]) == ['F']
        assert sorted([f.id for f in cls.rec.feature_after(15)]) == ['F', 'H']
        assert sorted([f.id for f in cls.rec.feature_after(15, nearest=True)]) == ['H']
        assert sorted([f.id for f in cls.rec.feature_after(15, strand=1)]) == ['H']
        assert sorted([f.id for f in cls.rec.feature_after(15, strand=-1)]) == ['F']
        assert sorted([f.id for f in cls.rec.feature_after(4, strand=-1)]) == []

    @classmethod
    def test_features_before(cls):
        assert sorted([f.id for f in cls.rec.feature_before(16)]) == ['B', 'G']
        assert sorted([f.id for f in cls.rec.feature_before(16, nearest=True)]) == ['G']
        assert sorted([f.id for f in cls.rec.feature_before(16, strand=1)]) == ['B']
        assert sorted([f.id for f in cls.rec.feature_before(16, strand=-1)]) == ['G']
        assert sorted([f.id for f in cls.rec.feature_before(20)]) == ['H']
        assert sorted([f.id for f in cls.rec.feature_before(20, strand=1)]) == ['H']
        assert sorted([f.id for f in cls.rec.feature_before(20, strand=-1)]) == []

    @classmethod
    def test_surrounding(cls):
        assert sorted([f.id for f in cls.rec.surrounding_features(16)]) == ['B', 'D', 'F', 'G']
        assert sorted([f.id for f in cls.rec.surrounding_features(16, strand=1)]) == ['B', 'D']
        assert sorted([f.id for f in cls.rec.surrounding_features(16, strand=-1)]) == ['F', 'G']
        assert sorted([f.id for f in cls.rec.surrounding_features(15)]) == ['B', 'F', 'H']
        assert sorted([f.id for f in cls.rec.surrounding_features(4)]) == ['A', 'B', 'F']
        assert sorted([f.id for f in cls.rec.surrounding_features(26)]) == ['C', 'G']
        assert sorted([f.id for f in cls.rec.surrounding_features(26, nearest=True)]) == ['G']
        assert sorted([f.id for f in cls.rec.surrounding_features(16, nearest=True)]) == ['G']
        assert sorted([f.id for f in cls.rec.surrounding_features(16, nearest=True, strand=1)]) == ['D']
        assert sorted([f.id for f in cls.rec.surrounding_features(16, nearest=True, strand=-1)]) == ['G']
        assert sorted([f.id for f in cls.rec.surrounding_features(11, nearest=True)]) == ['B', 'F']
        assert sorted([f.id for f in cls.rec.surrounding_features(11, nearest=True, strand=1)]) == ['B']
        assert sorted([f.id for f in cls.rec.surrounding_features(11, nearest=True, strand=-1)]) == ['F']

class SeqRecordStitchingTests(unittest.TestCase):
    def test_stitch(self):
        assert True is False


class SeqRecordComparisonTests(unittest.TestCase):
    rec1 = SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGC'), id='Rec1', name='DummyDNA')
    rec2 = SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGC'), id='Rec2', name='DummyDNA')
    rec3 = SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGC'), id='Rec3', name='DummyDNA')

    @classmethod
    def test_comparisons(cls):
        assert cls.rec1.__lt__(cls.rec2)
        assert not cls.rec3.__lt__(cls.rec2)
        assert cls.rec1.__le__(cls.rec2)
        assert cls.rec1.__le__(cls.rec1)
        assert not cls.rec2.__le__(cls.rec1)
        assert cls.rec1.__eq__(cls.rec1)
        assert not cls.rec3.__eq__(cls.rec2)
        assert cls.rec1.__ne__(cls.rec2)
        assert not cls.rec3.__ne__(cls.rec3)
        assert cls.rec2.__gt__(cls.rec1)
        assert not cls.rec1.__gt__(cls.rec2)
        assert cls.rec3.__ge__(cls.rec2)
        assert cls.rec1.__ge__(cls.rec1)
        assert not cls.rec1.__ge__(cls.rec2)
        assert cls.rec3 >= cls.rec3 > cls.rec1 == cls.rec1 < cls.rec2 <= cls.rec2
        assert sorted([cls.rec3, cls.rec1, cls.rec2]) == [cls.rec1, cls.rec2, cls.rec3]