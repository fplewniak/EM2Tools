"""
Tests for the extension module to the Biopython Bio.SeqUtils module.
"""
import unittest
import pytest
from Bio.SeqFeature import FeatureLocation

from EM2Libs import seq_utils
from EM2Libs.seq_feature import SeqFeatureEM2
from EM2Libs.seq_utils import SeqFilter as SF
from EM2Libs.seq_utils import GFF
from EM2Libs.seq import SeqEM2
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from pandas.testing import assert_frame_equal


class Pattern2Regex(unittest.TestCase):
    """Testing SeqUtils functions"""
    ambiguous_dna = {'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT',
                     'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
    ambiguous_prot = {'B': 'DN', 'Z': 'EQ', 'X': '.'}

    @staticmethod
    def test_ambiguous2string_nt():
        for k in Pattern2Regex.ambiguous_dna:
            assert seq_utils.ambiguous2string(k) == Pattern2Regex.ambiguous_dna[k]
        for r in 'ACGT':
            assert seq_utils.ambiguous2string(r, protein=True) == r
        with pytest.raises(ValueError, match=r'Invalid nucleotide .*'):
            seq_utils.ambiguous2string('O')

    @staticmethod
    def test_ambiguous2string_prot():
        for k in Pattern2Regex.ambiguous_prot:
            assert seq_utils.ambiguous2string(k, protein=True) == Pattern2Regex.ambiguous_prot[k]
        for r in 'ACDEFGHIKLMNPQRSTVWWY':
            assert seq_utils.ambiguous2string(r, protein=True) == r
        with pytest.raises(ValueError, match=r'Invalid amino-acid .*'):
            seq_utils.ambiguous2string('O', protein=True)

    @staticmethod
    def test_isambiguous_nt():
        for k in Pattern2Regex.ambiguous_dna:
            assert seq_utils.isambiguous(k) is True
        for k in 'ACGT([]){}0123456789,<>':
            assert seq_utils.isambiguous(k) is False

    @staticmethod
    def test_isambiguous_prot():
        for k in Pattern2Regex.ambiguous_prot:
            assert seq_utils.isambiguous(k, protein=True) is True
        for k in 'ACDEFGHIKLMNPQRSTVWY([]){}0123456789,<>':
            assert seq_utils.isambiguous(k, protein=True) is False

    @staticmethod
    def test_pattern2regex_nt():
        assert seq_utils.pattern2regex(
            'GARYSWKMBDHVNCT') == 'GA[AG][CT][CG][AT][GT][AC][CGT][AGT][ACT][ACG][ACGT]CT'

    @staticmethod
    def test_pattern2regex_prot():
        assert seq_utils.pattern2regex('FREDZXBAR', protein=True) == 'FRED[EQ].[DN]AR'

    @staticmethod
    def test_pattern2regex_any_none():
        assert seq_utils.pattern2regex('FRED[ANY]{ECPT}', protein=True) == 'FRED[ANY][^ECPT]'

    @staticmethod
    def test_pattern2regex_ambiguous():
        assert seq_utils.pattern2regex('FRED[AB]{ZR}', protein=True) == 'FRED[ADN][^EQR]'

    @staticmethod
    def test_pattern2regex_repeat():
        assert seq_utils.pattern2regex('FREDR(2)G(2,5)(KAI)(3,)[ILMV](4)',
                                       protein=True) == 'FREDR{2}G{2,5}(KAI){3,}[ILMV]{4}'

    @staticmethod
    def test_pattern2regex_start():
        assert seq_utils.pattern2regex('<FRED', protein=True) == '^FRED'

    @staticmethod
    def test_pattern2regex_end():
        assert seq_utils.pattern2regex('FRED>', protein=True) == 'FRED$'


def test_filter_by_length(dna_records):
    assert [r.id for r in SF().length(minlength=7).apply(dna_records)] == ['DNA1', 'DNA4', 'DNA5']
    assert [r.id for r in SF().length(7).apply(dna_records)] == ['DNA1', 'DNA4', 'DNA5']
    assert [r.id for r in SF().length(6, 10).apply(dna_records)] == ['DNA3', 'DNA4']
    assert SF().length(100).apply(dna_records) == []


def test_filter_by_pattern(dna_records):
    assert [r.id for r in SF().pattern('CATGW').apply(dna_records)] == ['DNA1', 'DNA5']
    assert SF().pattern('GGGGGGG').apply(dna_records) == []


def test_filter_by_name(dna_records):
    assert [r.id for r in SF().name(r'.*[35]$').apply(dna_records)] == ['DNA3', 'DNA5']
    assert SF().name(r'.*[6]$').apply(dna_records) == []


def test_filter_multi(dna_records):
    assert [r.id for r in SF().pattern('ACA').name(r'.*[1]$').length(7).apply(dna_records)] == \
           ['DNA1']
    assert SF().pattern('ACA').name(r'.*[2]$').length(7).apply(dna_records) == []


def test_by_length_keep_false(dna_records):
    assert [r.id for r in SF().length(7).keep(False).apply(dna_records)] == ['DNA2', 'DNA3']
    assert [r.id for r in SF().length(maxlength=6).apply(dna_records)] == ['DNA2', 'DNA3']


def test_no_filter(dna_records):
    assert SF().keep(False).apply(dna_records) == SF().apply(dna_records)


def test_by_pattern_keep_false(dna_records):
    assert [r.id for r in SF().pattern('CAGW').keep(False).apply(dna_records)] == ['DNA2']


def test_by_name_keep_false(dna_records):
    assert [r.id for r in SF().name(r'.*[345]$').keep(False).apply(dna_records)] == ['DNA1', 'DNA2']


class GFFtests(unittest.TestCase):
    sprot: SeqRecord = SeqRecord(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGC'), id='X',
                                 name='DummyDNA')
    sprot.features = [
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(0, 2), type='start', strand=1,
                      qualifiers={'codon': 'start', 'source': '', 'phase': '0', 'score': '0'}),
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(8, 18), type='domain', id='d1',
                      strand=0,
                      qualifiers={'source': '', 'phase': '0', 'score': '0'}),
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(16, 30), type='domain', id='d2',
                      strand=-1,
                      qualifiers={'source': '', 'phase': '0', 'score': '0'})
        ]

    df0 = DataFrame(
        {'seq_id': ['X'], 'source': [''], 'type': ['start'], 'start': ['0'], 'end': ['2'],
         'score': ['0'],
         'strand': ['+'], 'phase': ['0'], 'attributes': ['codon=start;id=<unknown id>']})

    df1 = DataFrame(
        {'seq_id': ['X', 'X', 'X'], 'source': ['', '', ''],
         'type': ['start', 'domain', 'domain'],
         'start': ['0', '8', '16'], 'end': ['2', '18', '30'], 'score': ['0', '0', '0'],
         'strand': ['+', '?', '-'], 'phase': ['0', '0', '0'],
         'attributes': ['codon=start;id=<unknown id>', 'id=d1', 'id=d2']})

    @classmethod
    def test_df_from_feature(cls):
        assert_frame_equal(GFF.df_from_feature(cls.sprot.features[0]).reset_index(drop=True),
                           cls.df0.reset_index(drop=True))
        assert_frame_equal(GFF.df_from_feature(None).reset_index(drop=True),
                           DataFrame(columns=['seq_id', 'source', 'type', 'start', 'end', 'score',
                                              'strand', 'phase',
                                              'attributes']).reset_index(drop=True))

    @classmethod
    def test_from_feature_list(cls):
        assert_frame_equal(GFF(cls.sprot.features[0:3]).df.reset_index(drop=True),
                           cls.df1.reset_index(drop=True))
        assert_frame_equal(GFF([]).df.reset_index(drop=True),
                           DataFrame(columns=['seq_id', 'source', 'type', 'start', 'end', 'score',
                                              'strand', 'phase',
                                              'attributes']).reset_index(drop=True))
        assert_frame_equal(GFF(None).df.reset_index(drop=True),
                           DataFrame(columns=['seq_id', 'source', 'type', 'start', 'end', 'score',
                                              'strand', 'phase',
                                              'attributes']).reset_index(drop=True))
        assert_frame_equal(GFF().df.reset_index(drop=True),
                           DataFrame(columns=['seq_id', 'source', 'type', 'start', 'end', 'score',
                                              'strand', 'phase',
                                              'attributes']).reset_index(drop=True))

    @classmethod
    def test_add_list(cls):
        assert_frame_equal(
            GFF([cls.sprot.features[0]]).add_feature_list(cls.sprot.features[1:3]).df.reset_index(
                drop=True),
            cls.df1.reset_index(drop=True))

    @classmethod
    def test_to_feature_list(cls):
        for i in range(0, len(cls.sprot.features)):
            assert GFF(input_df=cls.df1).to_feature_list(parents=cls.sprot)[i].__str__() == \
                   cls.sprot.features[i].__str__()
            assert GFF(input_df=cls.df1).to_feature_list(parents=cls.sprot)[i].parent.id == \
                   cls.sprot.features[i].parent.id
            assert GFF(input_df=cls.df1).to_feature_list()[i].parent is None
        for i in range(0, len(cls.sprot.features)):
            assert GFF(input_df=cls.df1).to_feature_list(parents=[cls.sprot,
                                                                  cls.sprot,
                                                                  cls.sprot])[i].__str__() == \
                   cls.sprot.features[i].__str__()
        with pytest.raises(ValueError,
                           match=r'The number of parents should match the number of features .*'):
            GFF(input_df=cls.df1).to_feature_list(parents=[cls.sprot])
