import unittest
from EM2Libs.Seq import SeqEM2


class SeqEM2Test(unittest.TestCase):
    @staticmethod
    def test_is_protein():
        assert SeqEM2.protein('HITHEREFRED').is_protein() == True
        assert SeqEM2.dna('CATAGTTAG').is_protein() == False


class SeqSearch(unittest.TestCase):
    @staticmethod
    def test_re_search():
        assert {m.group(): m.start() for m in SeqEM2.dna('CATAGCGCGCGGCTAAA').re_search(r'TAG')} == {'TAG': 2}
        assert {m.group(): m.start() for m in SeqEM2.dna('CATAGCGCGCGGCTAAA').re_search(r'CG{2}')} == {'CGG': 9}
        assert {m.group(): m.start() for m in SeqEM2.dna('CATAGCGCGCGGCTAAA').re_search(r'(GC){2,3}')} == {'GCGCGC': 4}
        assert {m.group(): m.start() for m in SeqEM2.dna('CATAGCGCGCGGCTAAA').re_search(r'[GC]{2,10}')} == {
            'GCGCGCGGC': 4}
        assert {m.group(): m.start() for m in SeqEM2.dna('CATAGCGCGCGGCTAAA').re_search(r'[^GC]{2,}')} == {'ATA': 1,
                                                                                                           'TAAA': 13}
        assert {m.group(): m.start() for m in SeqEM2.protein('HITHEREFREDANDGREG').re_search(r'THERE')} == {'THERE': 2}
        assert {m.group(): m.start() for m in SeqEM2.protein('HITHEREFREDANDGREG').re_search(r'(RE[FD]){2,3}')} == {
            'REFRED': 5}
        assert {m.group(): m.start() for m in SeqEM2.protein('HITHEREFREDANDGREG').re_search(r'[^GC]{2,}')} == {
            'HITHEREFREDAND': 0, 'RE': 15}

    @staticmethod
    def test_search():
        assert {m.group(): m.start() for m in SeqEM2.dna('CATAGCGCGCGGCTAAA').search('TAG')} == {'TAG': 2}
        assert {m.group(): m.start() for m in SeqEM2.dna('CATAGCGCGCGGCTAAA').search('CG(2)')} == {'CGG': 9}
        assert {m.group(): m.start() for m in SeqEM2.dna('CATAGCGCGCGGCTAAA').search('(GC)(2,3)')} == {'GCGCGC': 4}
        assert {m.group(): m.start() for m in SeqEM2.dna('CATAGCGCGCGGCTAAA').search('[GC](2,10)')} == {
            'GCGCGCGGC': 4}
        assert {m.group(): m.start() for m in SeqEM2.dna('CATAGCGCGCGGCTAAA').search('{GC}(2,)')} == {'ATA': 1,
                                                                                                      'TAAA': 13}
        assert {m.group(): m.start() for m in SeqEM2.protein('HITHEREFREDANDGREG').search('THERE')} == {'THERE': 2}
        assert {m.group(): m.start() for m in SeqEM2.protein('HITHEREFREDANDGREG').search('(RE[FD])(2,3)')} == {
            'REFRED': 5}
        assert {m.group(): m.start() for m in SeqEM2.protein('HITHEREFREDANDGREG').search('{GC}(2,)')} == {
            'HITHEREFREDAND': 0, 'RE': 15}
