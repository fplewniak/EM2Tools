#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#

import unittest

import pytest
from Bio.SeqFeature import FeatureLocation

from EM2Libs.Seq import SeqEM2
from EM2Libs.SeqFeature import SeqFeatureEM2
from EM2Libs.SeqRecord import SeqRecordEM2


class SeqRecord_FeatureTests(unittest.TestCase):
    rec: SeqRecordEM2 = SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGC'), id='Rec', name='DummyDNA')

    rec.features = [
        SeqFeatureEM2(parent=rec, location=FeatureLocation(0, 2), strand=1, id='A'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(6, 10), strand=1, id='B'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(15, 20), strand=1, id='C'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(20, 30), strand=1, id='D'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(6, 10), strand=-1, id='F'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(18, 25), strand=-1, id='G')
    ]

    @classmethod
    def test_overlap(cls):
        assert [f.id for f in cls.rec.overlap(16, 25)] == ['C', 'D', 'G']
        assert [f.id for f in cls.rec.overlap(16)] == ['C']

    def test_features_after(self):
        assert True is False

    def test_features_before(self):
        assert True is False

    def test_surrounding(self):
        assert True is False

class SeqRecord_StitchingTests(unittest.TestCase):
    def test_stitch(self):
        assert True is False
