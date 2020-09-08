#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
from em2lib.table import Table
from em2lib.table import TableTransform
from pandas import DataFrame
from pandas import Series
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
    assert TableTransform(table3)\
                .cond_transform(cond=lambda x: x<5 if isinstance(x, Number) else False, iftrue=f)\
                .cond_transform(cond=lambda x: x>2 if isinstance(x, Number) else False, iftrue=f, original=False)\
               .result().equals(table13)
    assert TableTransform(table3)\
                .cond_transform(cond=lambda x: x<5 if isinstance(x, Number) else False, iftrue=f)\
                .cond_transform(cond=lambda x: x>2 if isinstance(x, Number) else False, iftrue=f)\
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
    assert TableTransform(table2).combine(table4, lambda s1, s2: s1 if s1.sum() > s2.sum() else s2)\
        .result().equals(table4)
    assert TableTransform(table2).combine(table4, lambda s1, s2: s1 if s1.sum() < s2.sum() else s2)\
        .result().equals(table2)
    assert TableTransform(table2).combine(table4, lambda s1, s2: s1-s2, columns=['C', 'D'])\
        .result().equals(DataFrame(data={'A': ['a', 'a', 'y', 'z'], 'B': ['b', 'x', 'w', 'a'],
                                         'C': [-6, 0, 0, 0], 'D': [0, 0, 0, 0]}))
    assert TableTransform(table2[['C', 'D']]).combine(table4[['C', 'D']], f_df)\
        .result().equals(DataFrame(data={'C': {0: 'Good', 1: 'Bad', 2: 'Bad', 3: 'Bad'},
                                         'D': {0: 'Bad', 1: 'Bad', 2: 'Bad', 3: 'Bad'}}))

def test_normalize(table2):
    assert TableTransform(table2).normalize(columns=['C', 'D'], axis=1, norm='max')\
        .result().equals(DataFrame({'A': ['a', 'a', 'y', 'z'],
                                    'B': ['b', 'x', 'w', 'a'],
                                    'C': [0.125, 1., 1., 0.5],
                                    'D': [1., 1., 0.6, 1.]}))
    assert (TableTransform(table2).normalize(columns=['C', 'D'], axis=1, norm='l1')\
        .result()[['C', 'D']].to_numpy() == normalize([[1, 8], [2, 2], [5, 3], [3, 6]], axis=1, norm='l1')).all()
    assert (TableTransform(table2).normalize(columns=['C', 'D'], axis=1, norm='l2')\
        .result()[['C', 'D']].to_numpy() == normalize([[1, 8], [2, 2], [5, 3], [3, 6]], axis=1, norm='l2')).all()
    assert (TableTransform(table2).normalize(columns=['C', 'D'], axis=0, norm='l2')\
        .result()[['C', 'D']].to_numpy() == normalize([[1, 8], [2, 2], [5, 3], [3, 6]], axis=0, norm='l2')).all()
    assert TableTransform(table2).normalize(columns=['C', 'D'], norm='max')\
        .result().equals(DataFrame({'A': ['a', 'a', 'y', 'z'],
                                    'B': ['b', 'x', 'w', 'a'],
                                    'C': [0.125, 0.25, 0.625, 0.375],
                                    'D': [1., 0.25, 0.375, 0.75]}))
