from EM2Libs.Seq import SeqEM2
from EM2Libs.SeqRecord import SeqRecordEM2
from EM2Libs.SeqFeature import SeqFeatureEM2
from Bio.SeqFeature import FeatureLocation

import pytest

@pytest.fixture(scope="session")
def protein_sequence():
    return SeqEM2.protein('HITHEREFREDANDGREG')


@pytest.fixture(scope="session")
def dna_sequence():
    return SeqEM2.dna('CATAGCGCGCGGCTAAA')


@pytest.fixture(scope="session")
def protein_record():
    return SeqRecordEM2(SeqEM2.protein('HITHEREFREDANDGREG'))


@pytest.fixture(scope="session")
def dna_rec():
    rec: SeqRecordEM2 = SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGC'), id='Rec', name='DummyDNA')
    rec.features = [
        SeqFeatureEM2(parent=rec, location=FeatureLocation(0, 2), strand=1, id='A'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(6, 10), strand=1, id='B'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(15, 20), strand=1, id='C'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(20, 30), strand=1, id='D'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(6, 10), strand=-1, id='F'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(18, 25), strand=-1, id='G'),
        SeqFeatureEM2(parent=rec, location=FeatureLocation(16, 19), strand=0, id='H')
    ]
    return rec


@pytest.fixture(scope="session")
def dna_rec1():
    rec1 = SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGC'), id='Rec1', name='DummyDNA')
    rec1.features = [
        SeqFeatureEM2(parent=rec1, location=FeatureLocation(0, 5), strand=1, id='A1'),
        SeqFeatureEM2(parent=rec1, location=FeatureLocation(28, 33), strand=1, id='B1')
    ]
    return rec1

@pytest.fixture(scope="session")
def dna_rec2():
    rec2 = SeqRecordEM2(SeqEM2.dna('CAGCTGACGCATGAGTCGGTAACGATGCATGCATG'), id='Rec2', name='DummyDNA')
    rec2.features = [
        SeqFeatureEM2(parent=rec2, location=FeatureLocation(0, 5), strand=1, id='A2'),
        SeqFeatureEM2(parent=rec2, location=FeatureLocation(28, 33), strand=1, id='B2')
    ]
    return rec2

@pytest.fixture(scope="session")
def dna_rec3():
    rec3 = SeqRecordEM2(SeqEM2.dna('CACCTGACGCATGAGTCGGTAACGATGCATGCATG'), id='Rec3', name='DummyDNA')
    return rec3
