#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
from em2lib.table import Table


def test_table_read_csv(table1, table1_as_string, table1_as_dict):
    assert Table.from_csv("data/table1.csv", sep='\t').to_string() == table1_as_string
    assert Table.from_csv("data/table1.csv", sep='\t').to_dict() == table1_as_dict


def test_table_common_keys(table2, table3, table_2_3_common_keys):
    assert table2.get_common_keys(table3, index_self=['A', 'B'],
                                  index_other=['E', 'G']).to_string() == table_2_3_common_keys.to_string()


def test_table_different_keys(table2, table3, table_in_2_not_in_3, table_in_3_not_in_2):
    assert table2.get_keys_not_in(table3, index_self=['A', 'B'],
                                  index_other=['E', 'G']).to_string() == table_in_2_not_in_3.to_string()
    assert table3.get_keys_not_in(table2, index_self=['E', 'G'],
                                  index_other=['A', 'B']).to_string() == table_in_3_not_in_2.to_string()


def test_table_common_rows(table2, table3, table_2_2_common_rows, table_2_3_common_rows):
    assert table2.get_common_rows(table3, index_self=['A', 'B'], index_other=['E', 'G'],
                                  drop=True).to_string() == table_2_3_common_rows.to_string()
    assert table2.get_common_rows(table2, index_self=['A', 'B'], index_other=['A', 'B'], rsuffix='_2',
                                  drop=True).to_string() == table_2_2_common_rows.to_string()


def test_table_common_rows_full(table2, table3, table_2_3_common_rows_full):
    assert table2.get_common_rows(table3, index_self=['A', 'B'],
                                  index_other=['E', 'G'], drop=False).to_dict() == table_2_3_common_rows_full.to_dict()
