#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
from em2lib.table import Table
from em2lib.table import TableTransform
import em2lib.utils as em2
from pandas import DataFrame
from pandas import Series
from pandas import Index
from numbers import Number
import numpy as np
from sklearn.preprocessing import normalize


def test_table_common_keys(table2, table3, table_2_3_common_keys):
    assert Table.get_common_keys(table2, table3, keys_first=['A', 'B'],
                                 keys_other=['E', 'G']).to_string() == table_2_3_common_keys.to_string()


def test_table_different_keys(table2, table3, table_in_2_not_in_3, table_in_3_not_in_2):
    assert Table.get_keys_not_in(table2, table3, keys_first=['A', 'B'],
                                 keys_other=['E', 'G']).to_string() == table_in_2_not_in_3.to_string()
    assert Table.get_keys_not_in(table3, table2, keys_first=['E', 'G'],
                                 keys_other=['A', 'B']).to_string() == table_in_3_not_in_2.to_string()


def test_table_common_rows(table2, table3, table_2_2_common_rows, table_2_3_common_rows):
    assert Table.get_common_rows(table2, table3, keys_first=['A', 'B'], keys_other=['E', 'G'],
                                 drop=True).to_string() == table_2_3_common_rows.to_string()
    assert Table.get_common_rows(table2, table2, keys_first=['A', 'B'], keys_other=['A', 'B'], rsuffix='_2',
                                 drop=True).to_string() == table_2_2_common_rows.to_string()


def test_table_common_rows_full(table2, table3, table_2_3_common_rows_full):
    assert Table.get_common_rows(table2, table3, keys_first=['A', 'B'],
                                 keys_other=['E', 'G'], drop=False).to_dict() == table_2_3_common_rows_full.to_dict()


def test_table_rows_not_in(table2, table3, table_row_in_2_not_in_3, table_row_in_3_not_in_2):
    assert Table.get_rows_not_in(table2, table3, keys_first=['A', 'B'],
                                 keys_other=['E', 'G']).to_dict() == table_row_in_2_not_in_3.to_dict()
    assert Table.get_rows_not_in(table3, table2, keys_first=['E', 'G'],
                                 keys_other=['A', 'B']).to_dict() == table_row_in_3_not_in_2.to_dict()


def test_table_statistics(table2, table4):
    assert Table.statistics(table2, groupby='A', columns=['C', 'D'], func=[sum, np.mean]) \
        .equals(DataFrame({('sum', 'C'): {'a': 3, 'y': 5, 'z': 3},
                           ('sum', 'D'): {'a': 10, 'y': 3, 'z': 6},
                           ('mean', 'C'): {'a': 1.5, 'y': 5.0, 'z': 3.0},
                           ('mean', 'D'): {'a': 5.0, 'y': 3.0, 'z': 6.0}}))
    assert Table.statistics(table2, groupby='A', columns='C', func=[sum, np.mean]) \
        .equals(DataFrame({('sum', 'C'): {'a': 3, 'y': 5, 'z': 3},
                           ('mean', 'C'): {'a': 1.5, 'y': 5.0, 'z': 3.0}}))
    assert Table.statistics(table2, groupby='A', columns=['C', 'D'], func=[min, max]) \
        .equals(DataFrame({('min', 'C'): {'a': 1, 'y': 5, 'z': 3},
                           ('min', 'D'): {'a': 2, 'y': 3, 'z': 6},
                           ('max', 'C'): {'a': 2, 'y': 5, 'z': 3},
                           ('max', 'D'): {'a': 8, 'y': 3, 'z': 6}}))
    assert Table.statistics(table2.append(table4), groupby=['A', 'B'], columns=['C', 'D'], func=[sum, np.mean]) \
        .equals(DataFrame({('sum', 'C'): {('a', 'b'): 8, ('a', 'x'): 4, ('y', 'w'): 10, ('z', 'a'): 6},
                           ('sum', 'D'): {('a', 'b'): 16, ('a', 'x'): 4, ('y', 'w'): 6, ('z', 'a'): 12},
                           ('mean', 'C'): {('a', 'b'): 4.0, ('a', 'x'): 2.0, ('y', 'w'): 5.0, ('z', 'a'): 3.0},
                           ('mean', 'D'): {('a', 'b'): 8.0, ('a', 'x'): 2.0, ('y', 'w'): 3.0, ('z', 'a'): 6.0}}))
    assert Table.statistics(table2, columns=['C', 'D'], func=[sum, np.mean]) \
        .equals(DataFrame({'C': {'sum': 11.0, 'mean': 2.75}, 'D': {'sum': 19.0, 'mean': 4.75}}))
    assert Table.statistics(table2[['C', 'D']], func=[sum, np.mean]) \
        .equals(DataFrame({'C': {'sum': 11.0, 'mean': 2.75}, 'D': {'sum': 19.0, 'mean': 4.75}}))


def test_implode(table_org, table_expB, table_expBC, table_expBC_noidx,
                 table_expanded, table_expanded_imp, table_collapsed_all):
    assert Table.implode(table_expB).to_string() == table_org.to_string()
    assert Table.implode(table_expBC).to_string() == table_org.to_string()
    assert Table.implode(table_expBC_noidx, index='A').to_string() == table_org.to_string()
    assert Table.implode(table_expanded, index=0).to_string() == table_expanded_imp.to_string()
    assert Table.implode(table_collapsed_all).to_dict() == table_collapsed_all.to_dict()
    assert Table.implode(DataFrame({'V': {0: ('A', 'x', 'a'), 1: ('B', 'y', 'b'),
                                          2: ('A', 'z', 'c'), 3: ('B', 'w', 'd')}})).to_string() == \
           DataFrame({'V': {0: ('A', 'x', 'a'), 1: ('B', 'y', 'b'),
                            2: ('A', 'z', 'c'), 3: ('B', 'w', 'd')}}).to_string()


def test_collapse(table_expanded, table_collapsed_all, table_collapsed_1, table_expB, table_expBC):
    assert Table.collapse(table_expanded, groupby=0).to_dict() == table_collapsed_all.to_dict()
    assert Table.collapse(table_expanded, groupby=0, columns=1).to_string() == table_collapsed_1.to_string()
    assert Table.collapse(table_expB, groupby='A', columns='B', name='B').to_string() == DataFrame(
        {'B': {'a': ['x', 'y', 'z'], 'b': ['x', 'z'], 'c': 'y'}}, index=Index(['a', 'b', 'c'], name='A')).to_string()
    assert Table.collapse(table_expBC, groupby=['A', 'B'], name='C').to_string() == DataFrame(
        {'C': {('a', 'x'): ['x', 'y'], ('a', 'y'): ['x', 'y'], ('a', 'z'): ['x', 'y'],
               ('b', 'x'): 'x', ('b', 'z'): 'x', ('c', 'y'): ['y', 'z']}},
        index=Index([('a', 'x'), ('a', 'y'), ('a', 'z'), ('b', 'x'), ('b', 'z'), ('c', 'y')],
                    name=('A', 'B'))).to_string()


def test_expand(table_expanded, table_collapsed_all):
    assert Table.expand(table_collapsed_all).to_dict() == DataFrame({0: {0: 'A', 1: 'A', 2: 'B', 3: 'B'},
                                                                     1: {0: 'x', 1: 'z', 2: 'y', 3: 'w'},
                                                                     2: {0: 'a', 1: 'c', 2: 'b', 3: 'd'}}).to_dict()
    assert Table.expand(table_collapsed_all, columns=['col', 5]).to_dict() == DataFrame({
        0: {0: 'A', 1: 'A', 2: 'B', 3: 'B'},
        'col': {0: 'x', 1: 'z', 2: 'y', 3: 'w'},
        5: {0: 'a', 1: 'c', 2: 'b', 3: 'd'}}).to_dict()
    assert Table.expand(table_collapsed_all, columns='K').to_dict() == DataFrame({
        0: {0: 'A', 1: 'A', 2: 'B', 3: 'B'},
        'K': {0: 'x', 1: 'z', 2: 'y', 3: 'w'},
        'col_0': {0: 'a', 1: 'c', 2: 'b', 3: 'd'}}).to_dict()


def test_table_from_edges(df_edges):
    assert Table.table_from_edges(df_edges, no_link=0).equals(DataFrame({'a': {'a': 0, 'c': 3, 'b': 0, 'd': 0, 'e': 0},
                                                                         'c': {'a': 0, 'c': 0, 'b': 0, 'd': 0, 'e': 0},
                                                                         'b': {'a': 0, 'c': 0, 'b': 0, 'd': 0, 'e': 0},
                                                                         'd': {'a': 1, 'c': 0, 'b': 6, 'd': 0, 'e': 0},
                                                                         'e': {'a': 3, 'c': 0, 'b': 0, 'd': 0, 'e': 0}}
                                                                        ))
    assert Table.table_from_edges(df_edges, no_link=0, graph=False).equals(DataFrame({'d': {'a': 1, 'c': 0, 'b': 6},
                                                                                      'e': {'a': 3, 'c': 0, 'b': 0},
                                                                                      'a': {'a': 0, 'c': 3, 'b': 0}}))
    assert Table.table_from_edges(df_edges, no_link=0, directed=False).equals(DataFrame(
        {'a': {'a': 0, 'c': 3, 'b': 0, 'd': 1, 'e': 3}, 'c': {'a': 3, 'c': 0, 'b': 0, 'd': 0, 'e': 0},
         'b': {'a': 0, 'c': 0, 'b': 0, 'd': 6, 'e': 0}, 'd': {'a': 1, 'c': 0, 'b': 6, 'd': 0, 'e': 0},
         'e': {'a': 3, 'c': 0, 'b': 0, 'd': 0, 'e': 0}}))
    assert Table.table_from_edges(df_edges, no_link=0, link=1, directed=False).equals(DataFrame(
        {'a': {'a': 0, 'c': 1, 'b': 0, 'd': 1, 'e': 1}, 'c': {'a': 1, 'c': 0, 'b': 0, 'd': 0, 'e': 0},
         'b': {'a': 0, 'c': 0, 'b': 0, 'd': 1, 'e': 0}, 'd': {'a': 1, 'c': 0, 'b': 1, 'd': 0, 'e': 0},
         'e': {'a': 1, 'c': 0, 'b': 0, 'd': 0, 'e': 0}}))
    assert Table.table_from_edges(df_edges, no_link=0, link=1, directed=False).\
        equals(Table.table_from_edges(df_edges[[0,1]], no_link=0, link=1, directed=False))


def test_edges_from_table(df_table):
    assert Table.edges_from_table(df_table, no_link=0).to_dict() == DataFrame([['c', 'a', 3], ['d', 'a', 1],
                                                                               ['e', 'a', 3], ['a', 'c', 3],
                                                                               ['d', 'b', 6], ['a', 'd', 1],
                                                                               ['b', 'd', 6], ['a', 'e', 3]]).to_dict()
    assert Table.edges_from_table(df_table, no_link=0, directed=False).to_dict() == DataFrame([['c', 'a', 3],
                                                                                               ['d', 'a', 1],
                                                                                               ['e', 'a', 3],
                                                                                               ['d', 'b', 6]]).to_dict()


def gt3(x):
    if isinstance(x, Number):
        return x > 3
    return False


def dgt3(df):
    return df.applymap(lambda x: gt3(x))


def le3(x):
    if isinstance(x, Number):
        return x <= 3
    return False


def dle3(df):
    return df.applymap(lambda x: le3(x))


def mysqrt(x):
    if isinstance(x, Number):
        return np.sqrt(x)
    return x


def f(x):
    if isinstance(x, Number):
        return -x
    return x


def test_cond_transform(table2, table3, table6, table7, table8, table9, table10, table13, table14):
    assert TableTransform(table2) \
        .cond_transform(cond=lambda x: gt3(x), iftrue=lambda x: str(x) + '>3') \
        .cond_transform(cond=lambda x: le3(x), iftrue=lambda x: str(x) + '<=3') \
        .result().equals(table6)
    assert TableTransform(table2) \
        .cond_transform(cond=dgt3(table2), iftrue=lambda x: str(x) + '>3') \
        .cond_transform(cond=dle3(table2), iftrue=lambda x: str(x) + '<=3') \
        .result().equals(table6)
    assert TableTransform(table2) \
        .cond_transform(cond=dgt3(table2), iftrue=lambda x: str(x) + '>3', columns='C') \
        .cond_transform(cond=dle3(table2), iftrue=lambda x: str(x) + '<=3', columns='C') \
        .result().equals(table7)
    assert TableTransform(table2).cond_transform(cond=table2.eval('C + D > 5'), iftrue=mysqrt, columns='D') \
        .result().equals(table8)
    assert TableTransform(table2).cond_transform(cond=table2.eval('C + D > 5'), iftrue=table2.eval('D = C*2 +D')) \
        .result().equals(table9)
    assert TableTransform(table2).cond_transform(cond='C + D > 5', iftrue='D = C*2 +D') \
        .result().equals(table9)
    assert TableTransform(table2) \
        .cond_transform(cond=table2.eval('D + C> 5'), iftrue=table2.eval('C = D - C')) \
        .update() \
        .cond_transform(cond=lambda x: x < 0 if isinstance(x, Number) else False, iftrue=lambda x: 'X') \
        .result().equals(table10)
    assert TableTransform(table3) \
        .cond_transform(cond=lambda x: x < 5 if isinstance(x, Number) else False, iftrue=f) \
        .cond_transform(cond=lambda x: x > 2 if isinstance(x, Number) else False, iftrue=f, original=False) \
        .result().equals(table13)
    assert TableTransform(table3) \
        .cond_transform(cond=lambda x: x < 5 if isinstance(x, Number) else False, iftrue=f) \
        .cond_transform(cond=lambda x: x > 2 if isinstance(x, Number) else False, iftrue=f) \
        .result().equals(table14)


def test_replace(table2, table3, table11, table12):
    assert TableTransform(table3).replace(to_replace='a', value='A').result().equals(table11)
    assert TableTransform(table3).replace(to_replace='a', value='A',
                                          columns='E').result().equals(table12)


def f_df(s1, s2):
    d = []
    for x in s1.index:
        if s1[x] < s2[x]:
            d.append('Good')
        else:
            d.append('Bad')
    return Series(d)


def test_combine(table2, table4):
    assert TableTransform(table2).combine(table4, lambda s1, s2: s1 if s1.sum() > s2.sum() else s2) \
        .result().equals(table4)
    assert TableTransform(table2).combine(table4, lambda s1, s2: s1 if s1.sum() < s2.sum() else s2) \
        .result().equals(table2)
    assert TableTransform(table2).combine(table4, lambda s1, s2: s1 - s2, columns=['C', 'D']) \
        .result().equals(DataFrame(data={'A': ['a', 'a', 'y', 'z'], 'B': ['b', 'x', 'w', 'a'],
                                         'C': [-6, 0, 0, 0], 'D': [0, 0, 0, 0]}))
    assert TableTransform(table2[['C', 'D']]).combine(table4[['C', 'D']], f_df) \
        .result().equals(DataFrame(data={'C': {0: 'Good', 1: 'Bad', 2: 'Bad', 3: 'Bad'},
                                         'D': {0: 'Bad', 1: 'Bad', 2: 'Bad', 3: 'Bad'}}))


def test_normalize(table2):
    assert TableTransform(table2).normalize(columns=['C', 'D'], by='row', norm=em2.norm_max) \
        .result().equals(DataFrame({'A': ['a', 'a', 'y', 'z'],
                                    'B': ['b', 'x', 'w', 'a'],
                                    'C': [0.125, 1., 1., 0.5],
                                    'D': [1., 1., 0.6, 1.]}))
    assert TableTransform(table2).normalize(columns=['C', 'D'], by='col', norm=em2.norm_max) \
        .result().equals(DataFrame({'A': ['a', 'a', 'y', 'z'],
                                    'B': ['b', 'x', 'w', 'a'],
                                    'C': [0.2, 0.4, 1., 0.6],
                                    'D': [1., 0.25, 0.375, 0.75]}))
    assert TableTransform(table2).normalize(columns='C', by='col', norm=em2.norm_max) \
        .result().equals(DataFrame({'A': ['a', 'a', 'y', 'z'],
                                    'B': ['b', 'x', 'w', 'a'],
                                    'C': [0.2, 0.4, 1., 0.6],
                                    'D': [8, 2, 3, 6]}))
    assert (TableTransform(table2).normalize(columns=['C', 'D'], by='row', norm=em2.norm_sum_to_one) \
            .result()[['C', 'D']].to_numpy() == normalize([[1, 8], [2, 2], [5, 3], [3, 6]], axis=1, norm='l1')).all()
    assert TableTransform(table2).normalize(columns=['C', 'D'], norm=em2.norm_max) \
        .result().equals(DataFrame({'A': ['a', 'a', 'y', 'z'],
                                    'B': ['b', 'x', 'w', 'a'],
                                    'C': [0.125, 0.25, 0.625, 0.375],
                                    'D': [1., 0.25, 0.375, 0.75]}))
    assert TableTransform(table2).normalize(columns=['C', 'D'], by='row', norm=em2.norm_sum_to_one)\
        .result().equals(TableTransform(table2).normalize(by='row', norm=em2.norm_sum_to_one).result())
