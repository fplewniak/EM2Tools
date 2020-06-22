import unittest

import pytest
import operator

from em2lib.seq_feature import SeqFeatureEM2
from em2lib.seq_feature import FeatureFilter as FF
from em2lib.seq import SeqEM2
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
    def test_overlaps(cls):
        assert cls.sprot.features[2].overlaps(20, 25)
        assert cls.sprot.features[2].overlaps(20, 40)
        assert cls.sprot.features[2].overlaps(20)
        assert not cls.sprot.features[2].overlaps(35)
        assert not cls.sprot.features[2].overlaps(2, 5)

    @classmethod
    def test_overlaps_fuzzy(cls):
        with pytest.warns(UserWarning):
            cls.sprot.features[4].overlaps(35)
            cls.sprot.features[5].overlaps(3)

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
        assert cls.sprot.features[4].intersect(cls.sprot.features[7]).location == FeatureLocation(34, 37, ref='X')
        assert cls.sprot.features[2].intersect(cls.sprot.features[3]).location == cls.sprot.features[6].location
        assert cls.sprot.features[1].intersect(cls.sprot.features[3]).location == FeatureLocation(8, 18, ref='X')

    @classmethod
    def test_intersect_errors(cls):
        with pytest.raises(ValueError, match=r'Undetermined .*'):
            cls.sprot.features[0].intersect(SeqFeatureEM2(location=FeatureLocation(30, 37)))

    @classmethod
    def test_intersect_fuzzy(cls):
        with pytest.warns(UserWarning):
            cls.sprot.features[5].intersect(cls.sprot.features[0])

    @classmethod
    def test_move(cls):
        assert cls.sprot.features[0].move(5).location == FeatureLocation(5, 16)


def test_filter_length(dna_rec):
    assert FF().length(7, None).apply(dna_rec.features) == [dna_rec.features[i] for i in [3, 5]]


def test_filter_strand(dna_rec):
    assert FF().strand(1).apply(dna_rec.features) == dna_rec.features[0:4]
    assert FF().strand(-1).apply(dna_rec.features) == dna_rec.features[4:6]
    assert FF().strand(0).apply(dna_rec.features) == [dna_rec.features[6]]


def test_filter_location(dna_rec):
    assert FF().apply(dna_rec.features) == dna_rec.features
    assert FF().covers(24, 28).apply(dna_rec.features) == [dna_rec.features[3]]
    assert FF().overlaps(17, 24).apply(dna_rec.features) == [dna_rec.features[i] for i in [2,3,5,6]]
    assert FF().lies_within(12, 27).apply(dna_rec.features) == [dna_rec.features[i] for i in [2,5,6]]


def test_filter_type(dna_rec):
    assert FF().type('xxx').apply(dna_rec.features) == [dna_rec.features[5]]