from EM2Libs.Seq import SeqEM2
import pytest


@pytest.fixture(scope="class")
def seq_data(request):
    request.cls.prot = SeqEM2.protein('HITHEREFREDANDGREG')
    request.cls.dna = SeqEM2.dna('CATAGCGCGCGGCTAAA')

@pytest.fixture(scope="session")
def protein_sequence():
    return SeqEM2.protein('HITHEREFREDANDGREG')

@pytest.fixture(scope="session")
def dna_sequence():
    return SeqEM2.dna('CATAGCGCGCGGCTAAA')