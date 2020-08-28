#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
from em2lib.table import Table
from em2lib.table import TableTransform
from numbers import Number
import numpy as np


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
    else:
        return False


def dgt3(df):
    return df.applymap(lambda x: gt3(x))


def le3(x):
    if isinstance(x, Number):
        return x <= 3
    else:
        return False


def dle3(df):
    return df.applymap(lambda x: le3(x))


def mysqrt(x):
    if isinstance(x, Number):
        return np.sqrt(x)
    else:
        return x


def test_cond_transform(table2, table6, table7, table8, table9):
    assert TableTransform(table2).cond_transform(cond=lambda x: gt3(x), iftrue=lambda x: str(x) + '>3') \
               .cond_transform(cond=lambda x: le3(x),
                               iftrue=lambda x: str(x) + '<=3').result().to_string() == table6.to_string()
    assert TableTransform(table2).cond_transform(cond=dgt3(table2), iftrue=lambda x: str(x) + '>3') \
               .cond_transform(cond=dle3(table2),
                               iftrue=lambda x: str(x) + '<=3').result().to_string() == table6.to_string()
    assert TableTransform(table2).cond_transform(cond=dgt3(table2), iftrue=lambda x: str(x) + '>3', columns='C') \
               .cond_transform(cond=dle3(table2), iftrue=lambda x: str(x) + '<=3',
                               columns='C').result().to_string() == table7.to_string()
    assert TableTransform(table2).cond_transform(cond=table2.eval('C + D > 5'), iftrue=mysqrt,
                                                 columns='D').result().to_string() == table8.to_string()
    assert TableTransform(table2).cond_transform(cond=table2.eval('C + D > 5'), iftrue=table2.eval('D = C*2 +D')).\
               result().to_string() == table9.to_string()
