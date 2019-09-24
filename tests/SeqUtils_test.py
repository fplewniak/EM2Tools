"""
Tests for the extension module to the Biopython Bio.SeqUtils module.
"""
import unittest
import pytest
from EM2Libs import SeqUtils


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
        assert SeqUtils.pattern2regex('FREDR(2)G(2,5)(KAI)(3)[ILMV](4)',
                                      protein=True) == 'FREDR{2}G{2,5}(KAI){3}[ILMV]{4}'

    @staticmethod
    def test_pattern2regex_start():
        assert SeqUtils.pattern2regex('<FRED', protein=True) == '^FRED'

    @staticmethod
    def test_pattern2regex_end():
        assert SeqUtils.pattern2regex('FRED>', protein=True) == 'FRED$'
