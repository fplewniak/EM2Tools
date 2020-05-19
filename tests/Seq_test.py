import unittest
from EM2Libs.Seq import SeqEM2
import pytest


@pytest.fixture(scope="class")
def seq_data(request):
    request.cls.prot = SeqEM2.protein('HITHEREFREDANDGREG')
    request.cls.dna = SeqEM2.dna('CATAGCGCGCGGCTAAA')


@pytest.mark.usefixtures("seq_data")
class SeqEM2Test(unittest.TestCase):
    def test_is_protein(self):
        assert self.prot.is_protein() == True
        assert self.dna.is_protein() == False


@pytest.mark.usefixtures("seq_data")
class SeqSearch(unittest.TestCase):
    def test_re_search(self):
        assert {m.group(): m.start() for m in self.dna.re_search(r'TAG')} == {'TAG': 2}
        assert {m.group(): m.start() for m in self.dna.re_search(r'CG{2}')} == {'CGG': 9}
        assert {m.group(): m.start() for m in self.dna.re_search(r'(GC){2,3}')} == {'GCGCGC': 4}
        assert {m.group(): m.start() for m in self.dna.re_search(r'[GC]{2,10}')} == {'GCGCGCGGC': 4}
        assert {m.group(): m.start() for m in self.dna.re_search(r'[^GC]{2,}')} == {'ATA': 1, 'TAAA': 13}
        assert {m.group(): m.start() for m in self.prot.re_search(r'THERE')} == {'THERE': 2}
        assert {m.group(): m.start() for m in self.prot.re_search(r'(RE[FD]){2,3}')} == {'REFRED': 5}
        assert {m.group(): m.start() for m in self.prot.re_search(r'[^GC]{2,}')} == {'HITHEREFREDAND': 0, 'RE': 15}

    def test_search(self):
        assert {m.group(): m.start() for m in self.dna.search('TAG')} == {'TAG': 2}
        assert {m.group(): m.start() for m in self.dna.search('CG(2)')} == {'CGG': 9}
        assert {m.group(): m.start() for m in self.dna.search('(GC)(2,3)')} == {'GCGCGC': 4}
        assert {m.group(): m.start() for m in self.dna.search('[GC](2,10)')} == {'GCGCGCGGC': 4}
        assert {m.group(): m.start() for m in self.dna.search('{GC}(2,)')} == {'ATA': 1, 'TAAA': 13}
        assert {m.group(): m.start() for m in self.prot.search('THERE')} == {'THERE': 2}
        assert {m.group(): m.start() for m in self.prot.search('(RE[FD])(2,3)')} == {'REFRED': 5}
        assert {m.group(): m.start() for m in self.prot.search('{GC}(2,)')} == {'HITHEREFREDAND': 0, 'RE': 15}
        assert self.prot.search('RRRRRRR') == []
