"""
Tests for the extension module to the Biopython Bio.SeqUtils module.
"""
import unittest
import pytest
from Bio.SeqFeature import FeatureLocation

from EM2Libs import SeqUtils
from EM2Libs.SeqFeature import SeqFeatureEM2
from EM2Libs.SeqUtils import SeqFilter as SF
from EM2Libs.SeqUtils import GFF
from EM2Libs.Seq import SeqEM2
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from pandas.util.testing import assert_frame_equal


class Pattern2Regex(unittest.TestCase):
    """Testing SeqUtils functions"""
    ambiguous_dna = {'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
                     'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
    ambiguous_prot = {'B': 'DN', 'Z': 'EQ', 'X': '.'}

    @staticmethod
    def test_ambiguous2string_nt():
        for k in Pattern2Regex.ambiguous_dna:
            assert SeqUtils.ambiguous2string(k) == Pattern2Regex.ambiguous_dna[k]
        for r in 'ACGT':
            assert SeqUtils.ambiguous2string(r, protein=True) == r
        with pytest.raises(ValueError, match=r'Invalid nucleotide .*'):
            SeqUtils.ambiguous2string('O')

    @staticmethod
    def test_ambiguous2string_prot():
        for k in Pattern2Regex.ambiguous_prot:
            assert SeqUtils.ambiguous2string(k, protein=True) == Pattern2Regex.ambiguous_prot[k]
        for r in 'ACDEFGHIKLMNPQRSTVWWY':
            assert SeqUtils.ambiguous2string(r, protein=True) == r
        with pytest.raises(ValueError, match=r'Invalid amino-acid .*'):
            SeqUtils.ambiguous2string('O', protein=True)

    @staticmethod
    def test_isambiguous_nt():
        for k in Pattern2Regex.ambiguous_dna:
            assert SeqUtils.isambiguous(k) is True
        for k in 'ACGT([]){}0123456789,<>':
            assert SeqUtils.isambiguous(k) is False

    @staticmethod
    def test_isambiguous_prot():
        for k in Pattern2Regex.ambiguous_prot:
            assert SeqUtils.isambiguous(k, protein=True) is True
        for k in 'ACDEFGHIKLMNPQRSTVWY([]){}0123456789,<>':
            assert SeqUtils.isambiguous(k, protein=True) is False

    @staticmethod
    def test_pattern2regex_nt():
        assert SeqUtils.pattern2regex('GARYSWKMBDHVNCT') == 'GA[AG][CT][CG][AT][GT][AC][CGT][AGT][ACT][ACG][ACGT]CT'

    @staticmethod
    def test_pattern2regex_prot():
        assert SeqUtils.pattern2regex('FREDZXBAR', protein=True) == 'FRED[EQ].[DN]AR'

    @staticmethod
    def test_pattern2regex_any_none():
        assert SeqUtils.pattern2regex('FRED[ANY]{ECPT}', protein=True) == 'FRED[ANY][^ECPT]'

    @staticmethod
    def test_pattern2regex_ambiguous():
        assert SeqUtils.pattern2regex('FRED[AB]{ZR}', protein=True) == 'FRED[ADN][^EQR]'

    @staticmethod
    def test_pattern2regex_repeat():
        assert SeqUtils.pattern2regex('FREDR(2)G(2,5)(KAI)(3,)[ILMV](4)',
                                      protein=True) == 'FREDR{2}G{2,5}(KAI){3,}[ILMV]{4}'

    @staticmethod
    def test_pattern2regex_start():
        assert SeqUtils.pattern2regex('<FRED', protein=True) == '^FRED'

    @staticmethod
    def test_pattern2regex_end():
        assert SeqUtils.pattern2regex('FRED>', protein=True) == 'FRED$'


class SeqFilterTest(unittest.TestCase):
    records = [
        SeqRecord(SeqEM2.dna('ACAGTACCATGTAA'), id='DNA1', name='DNA1'),
        SeqRecord(SeqEM2.dna('ACAG'), id='DNA2', name='DNA2'),
        SeqRecord(SeqEM2.dna('ACAGTA'), id='DNA3', name='DNA3'),
        SeqRecord(SeqEM2.dna('ACAGTACCAT'), id='DNA4', name='DNA4'),
        SeqRecord(SeqEM2.dna('ACAGTACCATGT'), id='DNA5', name='DNA5')
    ]

    @classmethod
    def test_filter_by_length(cls):
        assert [r.id for r in SF().length(minlength=7).apply(cls.records)] == ['DNA1', 'DNA4', 'DNA5']
        assert [r.id for r in SF().length(7).apply(cls.records)] == ['DNA1', 'DNA4', 'DNA5']
        assert [r.id for r in SF().length(6, 10).apply(cls.records)] == ['DNA3', 'DNA4']
        assert SF().length(100).apply(cls.records) == []

    @classmethod
    def test_filter_by_pattern(cls):
        assert [r.id for r in SF().pattern('CATGW').apply(cls.records)] == ['DNA1', 'DNA5']
        assert SF().pattern('GGGGGGG').apply(cls.records) == []

    @classmethod
    def test_filter_by_name(cls):
        assert [r.id for r in SF().name(r'.*[35]$').apply(cls.records)] == ['DNA3', 'DNA5']
        assert SF().name(r'.*[6]$').apply(cls.records) == []

    @classmethod
    def test_filter_multi(cls):
        assert [r.id for r in SF().pattern('ACA').name(r'.*[1]$').length(7).apply(cls.records)] == ['DNA1']
        assert SF().pattern('ACA').name(r'.*[2]$').length(7).apply(cls.records) == []

    @classmethod
    def test_by_length_keep_false(cls):
        assert [r.id for r in SF().length(7).keep(False).apply(cls.records)] == ['DNA2', 'DNA3']
        assert [r.id for r in SF().length(maxlength=6).apply(cls.records)] == ['DNA2', 'DNA3']

    @classmethod
    def test_no_filter(cls):
        assert SF().keep(False).apply(cls.records) == SF().apply(cls.records)

    @classmethod
    def test_by_pattern_keep_false(cls):
        assert [r.id for r in SF().pattern('CAGW').keep(False).apply(cls.records)] == ['DNA2']

    @classmethod
    def test_by_name_keep_false(cls):
        assert [r.id for r in SF().name(r'.*[345]$').keep(False).apply(cls.records)] == ['DNA1', 'DNA2']


class GFFtests(unittest.TestCase):
    sprot: SeqRecord = SeqRecord(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGC'), id='X',
                                 name='DummyDNA')
    sprot.features = [
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(0, 2), type='start', strand=1,
                      qualifiers={'codon': 'start', 'source': '', 'phase': '0', 'score': '0'}),
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(8, 18), type='domain', id='d1', strand=0,
                      qualifiers={'source': '', 'phase': '0', 'score': '0'}),
        SeqFeatureEM2(parent=sprot, location=FeatureLocation(16, 30), type='domain', id='d2', strand=-1,
                      qualifiers={'source': '', 'phase': '0', 'score': '0'})
    ]

    df0 = DataFrame(
        {'seq_id': ['X'], 'source': [''], 'type': ['start'], 'start': ['0'], 'end': ['2'], 'score': ['0'],
         'strand': ['+'], 'phase': ['0'], 'attributes': ['codon=start;id=<unknown id>']})

    df1 = DataFrame(
        {'seq_id': ['X', 'X', 'X'], 'source': ['', '', ''], 'type': ['start', 'domain', 'domain'],
         'start': ['0', '8', '16'], 'end': ['2', '18', '30'], 'score': ['0', '0', '0'],
         'strand': ['+', '?', '-'], 'phase': ['0', '0', '0'],
         'attributes': ['codon=start;id=<unknown id>', 'id=d1', 'id=d2']})

    @classmethod
    def test_df_from_feature(cls):
        assert_frame_equal(GFF.df_from_feature(cls.sprot.features[0]).reset_index(drop=True),
                           cls.df0.reset_index(drop=True))
        assert_frame_equal(GFF.df_from_feature(None).reset_index(drop=True),
                           DataFrame(columns=['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                              'attributes']).reset_index(drop=True))

    @classmethod
    def test_from_feature_list(cls):
        assert_frame_equal(GFF(cls.sprot.features[0:3]).df.reset_index(drop=True), cls.df1.reset_index(drop=True))
        assert_frame_equal(GFF([]).df.reset_index(drop=True),
                           DataFrame(columns=['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                              'attributes']).reset_index(drop=True))
        assert_frame_equal(GFF(None).df.reset_index(drop=True),
                           DataFrame(columns=['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                              'attributes']).reset_index(drop=True))
        assert_frame_equal(GFF().df.reset_index(drop=True),
                           DataFrame(columns=['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                              'attributes']).reset_index(drop=True))

    @classmethod
    def test_addlist(cls):
        assert_frame_equal(
            GFF([cls.sprot.features[0]]).add_feature_list(cls.sprot.features[1:3]).df.reset_index(drop=True),
            cls.df1.reset_index(drop=True))

    @classmethod
    def test_to_feature_list(cls):
        for i in range(0, len(cls.sprot.features)):
            assert GFF(input_df=cls.df1).to_feature_list(parents=cls.sprot)[i].__str__() == cls.sprot.features[
                i].__str__()
        for i in range(0, len(cls.sprot.features)):
            assert GFF(input_df=cls.df1).to_feature_list(parents=[cls.sprot, cls.sprot, cls.sprot])[i].__str__() == \
                   cls.sprot.features[i].__str__()
        with pytest.raises(ValueError, match=r'The number of parents should match the number of features .*'):
            GFF(input_df=cls.df1).to_feature_list(parents=[cls.sprot])
