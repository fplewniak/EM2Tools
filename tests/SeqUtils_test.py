"""
Tests for the extension module to the Biopython Bio.SeqUtils module.
"""
import unittest
from EM2Libs import SeqUtils


class Pattern2Regex(unittest.TestCase):
    """Testing SeqUtils functions"""

    @staticmethod
    def test_ambiguous2string_nt():
        ambiguous_dna = {'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
                         'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
        for k in ambiguous_dna:
            assert SeqUtils.ambiguous2string(k) == ambiguous_dna[k]
        for r in 'ACGT':
            assert SeqUtils.ambiguous2string(r, protein=True) == r

    @staticmethod
    def test_ambiguous2string_prot():
        assert SeqUtils.ambiguous2string('B', protein=True) == 'DN'
        assert SeqUtils.ambiguous2string('Z', protein=True) == 'EQ'
        assert SeqUtils.ambiguous2string('X', protein=True) == '.'
        for r in 'ACDEFGHIKLMNPQRSTVWWY':
            assert SeqUtils.ambiguous2string(r, protein=True) == r

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
