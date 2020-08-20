#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
from em2lib.table import Table

def test_table_read_csv(table1, table1_as_string, table1_as_dict):
    assert Table.from_csv("data/table1.csv", sep='\t').to_string() == table1_as_string
    assert Table.from_csv("data/table1.csv", sep='\t').to_dict() == table1_as_dict

def test_table_common_keys(table2, table3, table_comon_keys):
    table2 = table2.set_index(['A', 'B'])
    table3 = table3.set_index(['E', 'G'])
    print(table2)
    assert table2.get_common_keys(table3).to_string() == table_comon_keys.to_string()

def test_table_different_keys(table2, table3, table_in_2_not_in_3, table_in_3_not_in_2):
    table2 = table2.set_index(['A', 'B'])
    table3 = table3.set_index(['E', 'G'])
    assert table2.get_keys_not_in(table3).to_string() == table_in_2_not_in_3.to_string()
    assert table3.get_keys_not_in(table2).to_string() == table_in_3_not_in_2.to_string()
