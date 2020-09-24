from em2lib.seq import SeqEM2
from em2lib.seq_record import SeqRecordEM2
from em2lib.seq_feature import SeqFeatureEM2
from Bio.SeqFeature import FeatureLocation

from pandas import DataFrame

from numpy import nan
from numpy import sqrt

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
# ##############################################Rec1_Rec2#######
def dna_seq(name='rec'):
    s = {'rec': 'ATGAGTCGGTAAAGATGCATGCGCGCCGCTGACGC',
         'rec1': 'ATGAGTCGGTAACGATGCATGCATGCAGCTGACGC',
         'rec2': 'CAGCTGACGCATGAGTCGGTAACGATGCATGCATG',
         'rec3': 'CACCTGACGCATGAGTCGGTAACGATGCATGCATG',
         'rec1NNNrec2': 'ATGAGTCGGTAACGATGCATGCATGCAGCTGACGCNNNCAGCTGACGCATGAGTCGGTAACGATGCATGCATG',
         'rec1_overlap_rec2/rec3': 'ATGAGTCGGTAACGATGCATGCATGCAGCTGACGCATGAGTCGGTAACGATGCATGCATG',
         'rec1_overlp_rec3_kpfls': 'ATGAGTCGGTAACGATGCATGCATGCACCTGACGCATGAGTCGGTAACGATGCATGCATG'
         }
    return SeqEM2.dna(s[name])


def joined_rec1NNNrec2():
    return SeqRecordEM2(dna_seq('rec1NNNrec2'), id='Rec1_NNN_Rec2', name='R1R2',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 5), strand=1, id='A1'),
                                  SeqFeatureEM2(location=FeatureLocation(28, 33), strand=1, id='B1'),
                                  SeqFeatureEM2(location=FeatureLocation(38, 43), strand=1, id='A2'),
                                  SeqFeatureEM2(location=FeatureLocation(66, 71), strand=1, id='B2')
                                  ])


def joined_rec1_overlap_rec2():
    return SeqRecordEM2(dna_seq('rec1_overlap_rec2/rec3'),
                        id='Rec1_overlap_Rec2', name='R1R2',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 5), strand=1, id='A1'),
                                  SeqFeatureEM2(location=FeatureLocation(28, 33), strand=1, id='B1'),
                                  SeqFeatureEM2(location=FeatureLocation(25, 30), strand=1, id='A2'),
                                  SeqFeatureEM2(location=FeatureLocation(53, 58), strand=1, id='B2')
                                  ])


@pytest.fixture(scope="session")
def dna_rec():
    return SeqRecordEM2(dna_seq('rec'), id='Rec', name='Rec',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 2), strand=1, id='A'),
                                  SeqFeatureEM2(location=FeatureLocation(6, 10), strand=1,
                                                type='yyy', id='B'),
                                  SeqFeatureEM2(location=FeatureLocation(15, 20), strand=1, id='C'),
                                  SeqFeatureEM2(location=FeatureLocation(20, 30), strand=1, id='D'),
                                  SeqFeatureEM2(location=FeatureLocation(6, 10), strand=-1, id='F'),
                                  SeqFeatureEM2(location=FeatureLocation(18, 25), strand=-1,
                                                type='xxx', id='G'),
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
def dna_rec2_rev():
    return SeqRecordEM2(dna_seq('rec2'), id='Rec2', name='R2',
                        features=[SeqFeatureEM2(location=FeatureLocation(0, 5), strand=1, id='A2'),
                                  SeqFeatureEM2(location=FeatureLocation(28, 33), strand=1, id='B2')
                                  ]).reverse_complement()


@pytest.fixture(scope="session")
def dna_rec3():
    return SeqRecordEM2(dna_seq('rec3'), id='Rec3', name='R3')


@pytest.fixture(scope="function")
def dna_rec1_NNN_rec2():
    return joined_rec1NNNrec2()


@pytest.fixture(scope="session")
def dna_stitch_rec1_NNN_rec2():
    stitched = joined_rec1NNNrec2()
    stitched.features += [SeqFeatureEM2(location=FeatureLocation(0, 35), strand=1, id='Rec1'),
                          SeqFeatureEM2(location=FeatureLocation(38, 73), strand=1, id='Rec2'),
                          SeqFeatureEM2(location=FeatureLocation(30, 48), strand=1, id='stitcher')
                          ]
    return stitched


@pytest.fixture(scope="session")
def dna_stitch_rec1_NNN_rec2rev():
    stitched = joined_rec1NNNrec2()
    stitched.features += [SeqFeatureEM2(location=FeatureLocation(0, 35), strand=1, id='Rec1'),
                          SeqFeatureEM2(location=FeatureLocation(38, 73), strand=-1, id='Rec2'),
                          SeqFeatureEM2(location=FeatureLocation(30, 48), strand=1, id='stitcher')
                          ]
    return stitched


@pytest.fixture(scope="session")
def dna_stitch_rec1_overlap_rec2():
    stitched = joined_rec1_overlap_rec2()
    stitched.features += [SeqFeatureEM2(location=FeatureLocation(0, 35), strand=1, id='Rec1'),
                          SeqFeatureEM2(location=FeatureLocation(25, 60), strand=1, id='Rec2'),
                          SeqFeatureEM2(location=FeatureLocation(22, 45), strand=1, id='stitcher')
                          ]
    return stitched


@pytest.fixture(scope="function")
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
        SeqRecordEM2(SeqEM2.dna('ACAGTACCATGTAA'), id='DNA1', name='DNA1'),
        SeqRecordEM2(SeqEM2.dna('ACAG'), id='DNA2', name='DNA2'),
        SeqRecordEM2(SeqEM2.dna('ACAGTA'), id='DNA3', name='DNA3'),
        SeqRecordEM2(SeqEM2.dna('ACAGTACCAT'), id='DNA4', name='DNA4'),
        SeqRecordEM2(SeqEM2.dna('ACAGTACCATGT'), id='DNA5', name='DNA5')
        ]


# #####################################################
# Examples of Tables and related data
# #####################################################
@pytest.fixture(scope="session")
def table1():
    return DataFrame(data={'A': {0: 4.0, 1: 6.2, 2: 1.3}, 'B': {0: 'x', 1: 'y', 2: 'z'}, 'C': {0: 8, 1: 9, 2: 5}})


@pytest.fixture(scope="session")
def table1_as_dict():
    return {'A': {0: 4.0, 1: 6.2, 2: 1.3}, 'B': {0: 'x', 1: 'y', 2: 'z'}, 'C': {0: 8, 1: 9, 2: 5}}


@pytest.fixture(scope="session")
def table1_as_string():
    return '     A  B  C\n0  4.0  x  8\n1  6.2  y  9\n2  1.3  z  5'


@pytest.fixture(scope="session")
def table2():
    return DataFrame(data={'A': ['a', 'a', 'y', 'z'], 'B': ['b', 'x', 'w', 'a'], 'C': [1, 2, 5, 3], 'D': [8, 2, 3, 6]})


@pytest.fixture(scope="session")
def table3():
    return DataFrame(data={'E': ['a', 'y', 'z', 'w'], 'F': [6, 8, 4, 2], 'G': ['b', 'w', 'b', 'a']})


@pytest.fixture(scope="session")
def table4():
    """
    :return: table with column C containing D - C value if D > C in table2
    """
    return DataFrame(data={'A': ['a', 'a', 'y', 'z'], 'B': ['b', 'x', 'w', 'a'], 'C': [7, 2, 5, 3], 'D': [8, 2, 3, 6]})


@pytest.fixture(scope="session")
def table5():
    """
    :return: table with column C containing D - C value if D > C in table2 and column D containing '*' if D + C > 8
    """
    return DataFrame(data={'A': ['a', 'a', 'y', 'z'], 'B': ['b', 'x', 'w', 'a'], 'C': [7, 2, 5, 3],
                           'D': ['X', 2, 3, 'X']})


@pytest.fixture(scope="session")
def table6():
    """
    :return: table2 copy with columns C and D modified to indicate whether the value is > 3 or <= 3
    """
    return DataFrame(data={'A': ['a', 'a', 'y', 'z'], 'B': ['b', 'x', 'w', 'a'], 'C': ['1<=3', '2<=3', '5>3', '3<=3'],
                           'D': ['8>3', '2<=3', '3<=3', '6>3']})


@pytest.fixture(scope="session")
def table7():
    """
    :return: table2 copy with columns C modified to indicate whether the value is > 3 or <= 3
    """
    return DataFrame(data={'A': ['a', 'a', 'y', 'z'], 'B': ['b', 'x', 'w', 'a'], 'C': ['1<=3', '2<=3', '5>3', '3<=3'],
                           'D': [8, 2, 3, 6]})


@pytest.fixture(scope="session")
def table8():
    """
     :return: table2 copy with columns D containing sqrt(x) if C + D > 5
     """
    return DataFrame(data={'A': ['a', 'a', 'y', 'z'], 'B': ['b', 'x', 'w', 'a'], 'C': [1, 2, 5, 3],
                           'D': [sqrt(8.0), 2.000000, sqrt(3.0), sqrt(6.0)]})


@pytest.fixture(scope="session")
def table9():
    """
     :return: table2 copy with columns D containing D = C*2 +D if C + D > 5
     """
    return DataFrame(data={'A': ['a', 'a', 'y', 'z'], 'B': ['b', 'x', 'w', 'a'], 'C': [1, 2, 5, 3],
                           'D': [10, 2, 13, 12]})


@pytest.fixture(scope="session")
def table10():
    """
    :return: table2 copy with columns C containing D - C if C + D > 5 or 'X' if the result was negative.
    """
    return DataFrame(
        data={'A': ['a', 'a', 'y', 'z'], 'B': ['b', 'x', 'w', 'a'], 'C': [7, 2, 'X', 3], 'D': [8, 2, 3, 6]})


@pytest.fixture(scope="session")
def table11():
    return DataFrame(data={'E': ['A', 'y', 'z', 'w'], 'F': [6, 8, 4, 2], 'G': ['b', 'w', 'b', 'A']})


@pytest.fixture(scope="session")
def table12():
    return DataFrame(data={'E': ['A', 'y', 'z', 'w'], 'F': [6, 8, 4, 2], 'G': ['b', 'w', 'b', 'a']})


@pytest.fixture(scope="session")
def table13():
    return DataFrame(data={'E': ['a', 'y', 'z', 'w'], 'F': [-6, -8, -4, -2], 'G': ['b', 'w', 'b', 'a']})


@pytest.fixture(scope="session")
def table14():
    return DataFrame(data={'E': ['a', 'y', 'z', 'w'], 'F': [-6, -8, 4, -2], 'G': ['b', 'w', 'b', 'a']})


table_original = DataFrame(data={'A': {0: 'a', 1: 'b', 2: 'c'},
                                 'B': {0: ['x', 'y', 'z'], 1: ['x', 'z'], 2: 'y'},
                                 'C': {0: ['x', 'y'], 1: 'x', 2: ['y', 'z']}})


@pytest.fixture(scope="session")
def table_org():
    return table_original


@pytest.fixture(scope="session")
def table_expB():
    return table_original.explode(column='B')


@pytest.fixture(scope="session")
def table_expBC():
    return table_original.explode(column='B').explode(column='C')


@pytest.fixture(scope="session")
def table_expBC_noidx():
    return DataFrame({'A': {0: 'a', 1: 'a', 2: 'a', 3: 'a', 4: 'a', 5: 'a', 6: 'b', 7: 'b', 8: 'c', 9: 'c'},
                      'B': {0: 'x', 1: 'x', 2: 'y', 3: 'y', 4: 'z', 5: 'z', 6: 'x', 7: 'z', 8: 'y', 9: 'y'},
                      'C': {0: 'x', 1: 'y', 2: 'x', 3: 'y', 4: 'x', 5: 'y', 6: 'x', 7: 'x', 8: 'y', 9: 'z'}})


@pytest.fixture(scope="session")
def table_expanded():
    return DataFrame({0: {0: 'A', 1: 'B', 2: 'A', 3: 'B'},
                      1: {0: 'x', 1: 'y', 2: 'z', 3: 'w'},
                      2: {0: 'a', 1: 'b', 2: 'c', 3: 'd'}})


@pytest.fixture(scope="session")
def table_expanded_imp():
    return DataFrame({0: {0: 'A', 1: 'B'}, 1: {0: ['x', 'z'], 1: ['y', 'w']}, 2: {0: ['a', 'c'], 1: ['b', 'd']}})


@pytest.fixture(scope="session")
def table_collapsed_all():
    t = DataFrame({tuple([1, 2]): {'A': [('x', 'a'), ('z','c')], 'B': [('y', 'b'), ('w','d')]}})
    t.index.rename(0, inplace=True)
    return t


@pytest.fixture(scope="session")
def table_collapsed_1():
    t = DataFrame({1: {'A': ['x', 'z'], 'B': ['y', 'w']}})
    t.index.rename(0, inplace=True)
    return t


@pytest.fixture(scope="session")
def df_edges():
    return DataFrame([['a', 'd', 1], ['a', 'e', 3],  ['c', 'a', 3], ['b', 'd', 6]])


@pytest.fixture(scope="session")
def df_table():
    return DataFrame({'a': {'a': 0, 'c': 3, 'b': 0, 'd': 1, 'e': 3}, 'c': {'a': 3, 'c': 0, 'b': 0, 'd': 0, 'e': 0},
                      'b': {'a': 0, 'c': 0, 'b': 0, 'd': 6, 'e': 0}, 'd': {'a': 1, 'c': 0, 'b': 6, 'd': 0, 'e': 0},
                      'e': {'a': 3, 'c': 0, 'b': 0, 'd': 0, 'e': 0}})


@pytest.fixture(scope="session")
def table_2_3_common_rows():
    return DataFrame(data={'A': ['a', 'y'], 'B': ['b', 'w'], 'C': [1, 5], 'D': [8, 3],
                           'F': [6, 8]})


@pytest.fixture(scope="session")
def table_2_3_common_rows_full():
    return DataFrame(data={'A': ['a', 'y'], 'B': ['b', 'w'], 'C': [1, 5], 'D': [8, 3],
                           'E': ['a', 'y'], 'F': [6, 8], 'G': ['b', 'w']})


@pytest.fixture(scope="session")
def table_2_2_common_rows():
    return DataFrame(data={'A': ['a', 'a', 'y', 'z'], 'B': ['b', 'x', 'w', 'a'], 'C': [1, 2, 5, 3], 'D': [8, 2, 3, 6],
                           'C_2': [1, 2, 5, 3], 'D_2': [8, 2, 3, 6]})


@pytest.fixture(scope="session")
def table_2_3_combined():
    return DataFrame(data={'A': ['a', 'a', 'y', 'z', nan, nan], 'B': ['b', 'x', 'w', 'a', nan, nan],
                           'C': [1, 2, 5, 3, nan, nan], 'D': [8, 2, 3, 6, nan, nan],
                           'E': ['a', nan, 'y', nan, 'z', 'w'], 'F': [6, nan, 8, nan, 4, 2],
                           'G': ['b', nan, 'w', nan, 'b', 'a']})


@pytest.fixture(scope="session")
def table_2_3_common_keys():
    return DataFrame(data={'A': ['a', 'y'], 'B': ['b', 'w']})


@pytest.fixture(scope="session")
def table_in_2_not_in_3():
    return DataFrame(data={'A': ['a', 'z'], 'B': ['x', 'a']})


@pytest.fixture(scope="session")
def table_in_3_not_in_2():
    return DataFrame(data={'E': ['z', 'w'], 'G': ['b', 'a']})


@pytest.fixture(scope="session")
def table_row_in_2_not_in_3():
    return DataFrame(data={'A': ['a', 'z'], 'B': ['x', 'a'], 'C': [2, 3], 'D': [2, 6]})


@pytest.fixture(scope="session")
def table_row_in_3_not_in_2():
    return DataFrame(data={'E': ['z', 'w'], 'F': [4, 2], 'G': ['b', 'a']})
