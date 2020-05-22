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
    return SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGC'), id='Rec', name='DummyDNA',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 2), strand=1, id='A'),
                                  SeqFeatureEM2(location=FeatureLocation(6, 10), strand=1, id='B'),
                                  SeqFeatureEM2(location=FeatureLocation(15, 20), strand=1, id='C'),
                                  SeqFeatureEM2(location=FeatureLocation(20, 30), strand=1, id='D'),
                                  SeqFeatureEM2(location=FeatureLocation(6, 10), strand=-1, id='F'),
                                  SeqFeatureEM2(location=FeatureLocation(18, 25), strand=-1, id='G'),
                                  SeqFeatureEM2(location=FeatureLocation(16, 19), strand=0, id='H')])


@pytest.fixture(scope="session")
def dna_rec1():
    return SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGC'), id='Rec1', name='DummyDNA',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 5), strand=1, id='A1'),
                                  SeqFeatureEM2(location=FeatureLocation(28, 33), strand=1, id='B1')])


@pytest.fixture(scope="session")
def dna_rec2():
    return SeqRecordEM2(SeqEM2.dna('CAGCTGACGCATGAGTCGGTAACGATGCATGCATG'), id='Rec2', name='DummyDNA',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 5), strand=1, id='A2'),
                                  SeqFeatureEM2(location=FeatureLocation(28, 33), strand=1, id='B2')])


@pytest.fixture(scope="session")
def dna_rec3():
    return SeqRecordEM2(SeqEM2.dna('CACCTGACGCATGAGTCGGTAACGATGCATGCATG'), id='Rec3', name='DummyDNA')
