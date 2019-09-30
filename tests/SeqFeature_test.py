import unittest

import pytest

from EM2Libs.SeqFeature import SeqFeatureEM2
from EM2Libs.Seq import SeqEM2
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import AfterPosition, BeforePosition


class SeqFeatureTests(unittest.TestCase):
    sprot: SeqRecord = SeqRecord(SeqEM2.protein('MYNAMEISFREDHEREIAMWHEREARETHEYALLTHISISEXCELLENT'), id='X',
                                 name='DummyProt')
    sprot.features = [
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(0, 11), type='domain', id='d1'),  # MYNAMEISFRED
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(8, 18), type='domain', id='d2'),  # FREDHEREIAM
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(19, 30), type='domain', id='d3'),  # WHEREARETHEY
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(6, 23), type='domain', id='d4'),  # ISFREDHEREIAMWHERE
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(34, AfterPosition(39)), id='d5'),  # THISIS
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(BeforePosition(2), 5), type='domain', id='d6'),  # MYNAME
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(19, 23), type='domain', id='d7'),  # WHERE
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(BeforePosition(30), 37), type='domain', id='d8')  # YALLTHI
    ]

    @classmethod
    def test_parent(cls):
        assert [f.id for f in cls.sprot.features] == ['d1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8']
        assert cls.sprot.features[1].parent.id == cls.sprot.id
        assert cls.sprot.features[1].parent.name == cls.sprot.name
        assert cls.sprot.features[1].parent.seq._data == cls.sprot.seq._data

    @classmethod
    def test_lies_within(cls):
        assert cls.sprot.features[1].lies_within(5, 25)
        assert not cls.sprot.features[1].lies_within(10, 25)
        assert not cls.sprot.features[1].lies_within(19, 25)

    @classmethod
    def test_lies_within_fuzzy(cls):
        with pytest.warns(UserWarning):
            cls.sprot.features[4].lies_within(30, 42)
            cls.sprot.features[5].lies_within(0, 10)


    @classmethod
    def test_contains(cls):
        assert cls.sprot.features[2].contains(20)
        assert not cls.sprot.features[2].contains(35)

    @classmethod
    def test_overlaps(cls):
        assert cls.sprot.features[2].overlaps(20, 25)
        assert cls.sprot.features[2].overlaps(20, 40)
        assert not cls.sprot.features[2].overlaps(2, 5)

    @classmethod
    def test_contains_fuzzy(cls):
        with pytest.warns(UserWarning):
            cls.sprot.features[4].contains(35)
            cls.sprot.features[5].contains(3)

    @classmethod
    def test_covers(cls):
        assert cls.sprot.features[3].covers(15, 20)
        assert not cls.sprot.features[3].covers(4, 20)

    @classmethod
    def test_covers_fuzzy(cls):
        with pytest.warns(UserWarning):
            cls.sprot.features[4].covers(35, 38)
            cls.sprot.features[5].covers(3, 4)

    @classmethod
    def test_intersect(cls):
        assert cls.sprot.features[4].intersect(cls.sprot.features[7]).location == FeatureLocation(34, 37)
        assert cls.sprot.features[2].intersect(cls.sprot.features[3]).location == cls.sprot.features[6].location
        assert cls.sprot.features[1].intersect(cls.sprot.features[3]).location == FeatureLocation(8, 18)

    @classmethod
    def test_intersect_errors(cls):
        with pytest.raises(ValueError, match=r'Undetermined .*'):
            cls.sprot.features[0].intersect(SeqFeatureEM2(location=FeatureLocation(30, 37)))

    @classmethod
    def test_intersect_fuzzy(cls):
        with pytest.warns(UserWarning):
            cls.sprot.features[5].intersect(cls.sprot.features[0])
