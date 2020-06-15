from em2lib.seq import SeqEM2
from em2lib.seq_record import SeqRecordEM2
from em2lib.seq_feature import SeqFeatureEM2
from Bio.SeqFeature import FeatureLocation

import pytest


# #####################################################
# Examples of protein sequences and records
# #####################################################
def protein_seq():
    return SeqEM2.protein('HITHEREFREDANDGREG')


@pytest.fixture(scope="session")
def protein_rec():
    return SeqRecordEM2(protein_seq())


# #####################################################

# #####################################################
# Examples of DNA sequences and records
# #####################################################
def dna_seq(name='rec1'):
    s = {'rec1': 'ATGAGTCGGTAACGATGCATGCATGCAGCTGACGC',
         'rec2': 'CAGCTGACGCATGAGTCGGTAACGATGCATGCATG',
         'rec3': 'CACCTGACGCATGAGTCGGTAACGATGCATGCATG',
         'rec1NNNrec2': 'ATGAGTCGGTAACGATGCATGCATGCAGCTGACGCNNNCAGCTGACGCATGAGTCGGTAACGATGCATGCATG',
         'rec1_overlap_rec2/rec3': 'ATGAGTCGGTAACGATGCATGCATGCAGCTGACGCATGAGTCGGTAACGATGCATGCATG',
         'rec1_overlp_rec3_kpfls': 'ATGAGTCGGTAACGATGCATGCATGCACCTGACGCATGAGTCGGTAACGATGCATGCATG'
         }
    return SeqEM2.dna(s[name])


def joined_rec1NNNrec2():
    return SeqRecordEM2(dna_seq('rec1NNNrec2'), id='Rec1_NNN_Rec2', name='R1R2',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 5), strand=1, id='A1',
                                                ref='<unknown id>'),
                                  SeqFeatureEM2(location=FeatureLocation(28, 33), strand=1, id='B1',
                                                ref='<unknown id>'),
                                  SeqFeatureEM2(location=FeatureLocation(38, 43), strand=1, id='A2',
                                                ref='<unknown id>'),
                                  SeqFeatureEM2(location=FeatureLocation(66, 71), strand=1, id='B2',
                                                ref='<unknown id>')
                                  ])


def joined_rec1_overlap_rec2():
    return SeqRecordEM2(dna_seq('rec1_overlap_rec2/rec3'),
                        id='Rec1_overlap_Rec2', name='R1R2',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 5), strand=1, id='A1',
                                                ref='<unknown id>'),
                                  SeqFeatureEM2(location=FeatureLocation(28, 33), strand=1, id='B1',
                                                ref='<unknown id>'),
                                  SeqFeatureEM2(location=FeatureLocation(25, 30), strand=1, id='A2',
                                                ref='<unknown id>'),
                                  SeqFeatureEM2(location=FeatureLocation(53, 58), strand=1, id='B2',
                                                ref='<unknown id>')
                                  ])


@pytest.fixture(scope="session")
def dna_rec():
    return SeqRecordEM2(dna_seq('rec1'), id='Rec', name='Rec',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 2), strand=1, id='A'),
                                  SeqFeatureEM2(location=FeatureLocation(6, 10), strand=1, id='B'),
                                  SeqFeatureEM2(location=FeatureLocation(15, 20), strand=1, id='C'),
                                  SeqFeatureEM2(location=FeatureLocation(20, 30), strand=1, id='D'),
                                  SeqFeatureEM2(location=FeatureLocation(6, 10), strand=-1, id='F'),
                                  SeqFeatureEM2(location=FeatureLocation(18, 25), strand=-1,
                                                id='G'),
                                  SeqFeatureEM2(location=FeatureLocation(16, 19), strand=0, id='H')
                                  ])


@pytest.fixture(scope="session")
def dna_rec1():
    return SeqRecordEM2(dna_seq('rec1'), id='Rec1', name='R1',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 5), strand=1, id='A1'),
                                  SeqFeatureEM2(location=FeatureLocation(28, 33), strand=1, id='B1')
                                  ])


@pytest.fixture(scope="session")
def dna_rec2():
    return SeqRecordEM2(dna_seq('rec2'), id='Rec2', name='R2',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 5), strand=1, id='A2'),
                                  SeqFeatureEM2(location=FeatureLocation(28, 33), strand=1, id='B2')
                                  ])


@pytest.fixture(scope="session")
def dna_rec3():
    return SeqRecordEM2(dna_seq('rec3'), id='Rec3', name='R3')


@pytest.fixture(scope="function")
def dna_rec1_NNN_rec2():
    return joined_rec1NNNrec2()


@pytest.fixture(scope="session")
def dna_stitch_rec1_NNN_rec2():
    stitched = joined_rec1NNNrec2()
    stitched.features += [SeqFeatureEM2(location=FeatureLocation(0, 34), strand=1, id='Rec1'),
                          SeqFeatureEM2(location=FeatureLocation(38, 73), strand=1, id='Rec2'),
                          SeqFeatureEM2(location=FeatureLocation(30, 48), strand=1, id='stitcher')
                          ]
    return stitched


@pytest.fixture(scope="session")
def dna_stitch_rec1_overlap_rec2():
    stitched = joined_rec1_overlap_rec2()
    stitched.features += [SeqFeatureEM2(location=FeatureLocation(0, 34), strand=1, id='Rec1'),
                          SeqFeatureEM2(location=FeatureLocation(25, 60), strand=1, id='Rec2'),
                          SeqFeatureEM2(location=FeatureLocation(22, 45), strand=1, id='stitcher')
                          ]
    return stitched


@pytest.fixture(scope="session")
def dna_rec1_overlap_rec2():
    return joined_rec1_overlap_rec2()


@pytest.fixture(scope="session")
def dna_rec1_overlap_rec3():
    return SeqRecordEM2(dna_seq('rec1_overlap_rec2/rec3'), id='Rec1_overlap_Rec3', name='R1R3')


@pytest.fixture(scope="session")
def dna_rec1_overlap_rec3_keep_false():
    return SeqRecordEM2(dna_seq('rec1_overlp_rec3_kpfls'), id='Rec1_overlp_Rec3_kpfls', name='R1R3')


@pytest.fixture(scope="session")
def dna_records():
    return [
        SeqRecordEM2(dna_seq('rec1'), id='DNA1', name='DNA1'),
        SeqRecordEM2(dna_seq('rec1'), id='DNA2', name='DNA2'),
        SeqRecordEM2(dna_seq('rec1'), id='DNA3', name='DNA3'),
        SeqRecordEM2(dna_seq('rec1'), id='DNA4', name='DNA4'),
        SeqRecordEM2(dna_seq('rec1'), id='DNA5', name='DNA5')
        ]
