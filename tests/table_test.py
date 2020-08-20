#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
from em2lib.table import Table

def test_table_read_csv(table1, table1_as_string, table1_as_dict):
    assert Table.from_csv("data/table1.csv", sep='\t').to_string() == table1_as_string
    assert Table.from_csv("data/table1.csv", sep='\t').to_dict() == table1_as_dict

