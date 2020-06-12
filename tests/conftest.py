from em2lib.seq import SeqEM2
from em2lib.seq_record import SeqRecordEM2
from em2lib.seq_feature import SeqFeatureEM2
from Bio.SeqFeature import FeatureLocation

import pytest


@pytest.fixture(scope="session")
def protein_rec():
    return SeqRecordEM2(SeqEM2.protein('HITHEREFREDANDGREG'))


@pytest.fixture(scope="session")
def dna_rec():
    return SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAAAGATGCATGCGCGCCGCTGACGC'), id='Rec', name='DummyDNA',
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
    return SeqRecordEM2(SeqEM2.dna('CACCTGACGCATGAGTCGGTAACGATGCATGCATG'), id='Rec3',
                        name='DummyDNA')


@pytest.fixture(scope="session")
def dna_rec1_NNN_rec2():
    return SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGCNNNCAGCTGACGCATGAGTCGGTAACGATGCATGCATG'),
                        id='Rec1_NNN_Rec2', name='DummyDNA',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 5), strand=1, id='A1', ref='<unknown id>'),
                                  SeqFeatureEM2(location=FeatureLocation(28, 33), strand=1, id='B1', ref='<unknown id>'),
                                  SeqFeatureEM2(location=FeatureLocation(38, 43), strand=1, id='A2', ref='<unknown id>'),
                                  SeqFeatureEM2(location=FeatureLocation(66, 71), strand=1, id='B2', ref='<unknown id>')
                                  ])


@pytest.fixture(scope="session")
def dna_rec1_overlap_rec2():
    return SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGCATGAGTCGGTAACGATGCATGCATG'),
                        id='Rec1_overlap_Rec2', name='DummyDNA',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 5), strand=1, id='A1', ref='<unknown id>'),
                                  SeqFeatureEM2(location=FeatureLocation(28, 33), strand=1, id='B1', ref='<unknown id>'),
                                  SeqFeatureEM2(location=FeatureLocation(25, 30), strand=1, id='A2', ref='<unknown id>'),
                                  SeqFeatureEM2(location=FeatureLocation(53, 58), strand=1, id='B2', ref='<unknown id>')
                                  ])


@pytest.fixture(scope="session")
def dna_rec1_overlap_rec3():
    return SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCAGCTGACGCATGAGTCGGTAACGATGCATGCATG'),
                        id='Rec1_overlap_Rec3', name='DummyDNA')


@pytest.fixture(scope="session")
def dna_rec1_overlap_rec3_keep_false():
    return SeqRecordEM2(SeqEM2.dna('ATGAGTCGGTAACGATGCATGCATGCACCTGACGCATGAGTCGGTAACGATGCATGCATG'),
                        id='Rec1_overlap_Rec3_keep_false', name='DummyDNA')


@pytest.fixture(scope="session")
def dna_records():
    return [
        SeqRecordEM2(SeqEM2.dna('ACAGTACCATGTAA'), id='DNA1', name='DNA1'),
        SeqRecordEM2(SeqEM2.dna('ACAG'), id='DNA2', name='DNA2'),
        SeqRecordEM2(SeqEM2.dna('ACAGTA'), id='DNA3', name='DNA3'),
        SeqRecordEM2(SeqEM2.dna('ACAGTACCAT'), id='DNA4', name='DNA4'),
        SeqRecordEM2(SeqEM2.dna('ACAGTACCATGT'), id='DNA5', name='DNA5')
    ]
