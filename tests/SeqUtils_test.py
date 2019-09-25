"""
Tests for the extension module to the Biopython Bio.SeqUtils module.
"""
import unittest
import pytest
from EM2Libs import SeqUtils
from EM2Libs.Seq import SeqEM2
from Bio.SeqRecord import SeqRecord


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


class SeqFilter(unittest.TestCase):
    records = [
        SeqRecord(SeqEM2.dna('ACAGTACCATGTAA'), id='DNA1', name='DNA1'),
        SeqRecord(SeqEM2.dna('ACAG'), id='DNA2', name='DNA2'),
        SeqRecord(SeqEM2.dna('ACAGTA'), id='DNA3', name='DNA3'),
        SeqRecord(SeqEM2.dna('ACAGTACCAT'), id='DNA4', name='DNA4'),
        SeqRecord(SeqEM2.dna('ACAGTACCATGT'), id='DNA5', name='DNA5')
    ]

    @staticmethod
    def test_filter_by_length():
        assert [r.id for r in SeqUtils.seqfilter(SeqFilter.records, minlength=7)] == ['DNA1', 'DNA4', 'DNA5']
        assert [r.id for r in SeqUtils.seqfilter(SeqFilter.records, minlength=6, maxlength=10)] == ['DNA3', 'DNA4']
        assert SeqUtils.seqfilter(SeqFilter.records, minlength=100) == []

    @staticmethod
    def test_filter_by_pattern():
        assert [r.id for r in SeqUtils.seqfilter(SeqFilter.records, pattern='CATGW')] == ['DNA1', 'DNA5']
        assert SeqUtils.seqfilter(SeqFilter.records, pattern='GGGGGGG') == []

    @staticmethod
    def test_filter_by_name():
        assert [r.id for r in SeqUtils.seqfilter(SeqFilter.records, name=r'.*[35]$')] == ['DNA3', 'DNA5']
        assert SeqUtils.seqfilter(SeqFilter.records, name=r'.*[6]$') == []

    @staticmethod
    def test_filter_multi():
        assert [r.id for r in SeqUtils.seqfilter(SeqFilter.records, pattern='ACA', name=r'.*[1]$', minlength=7)] == [
            'DNA1']
        assert SeqUtils.seqfilter(SeqFilter.records, pattern='ACA', name=r'.*[2]$', minlength=7) == []

    @staticmethod
    def test_by_length_keep_false():
        assert [r.id for r in SeqUtils.seqfilter(SeqFilter.records, minlength=7, keep=False)] == ['DNA2', 'DNA3']

    @staticmethod
    def test_no_filter():
        assert SeqUtils.seqfilter(SeqFilter.records, keep=False) == SeqUtils.seqfilter(SeqFilter.records)

    @staticmethod
    def test_by_pattern_keep_false():
        assert [r.id for r in SeqUtils.seqfilter(SeqFilter.records, pattern='CAGTW', keep=False)] == ['DNA2']

    @staticmethod
    def test_by_name_keep_false():
        assert [r.id for r in SeqUtils.seqfilter(SeqFilter.records, name=r'.*[345]$', keep=False)] == ['DNA1', 'DNA2']
