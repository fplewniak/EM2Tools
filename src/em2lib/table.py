"""
Extension of pandas.DataFrame class to add or make functionalities easier to use
"""
#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
from pandas import read_csv
from pandas import DataFrame


class Table(DataFrame):
    """
    Extension of pandas.DataFrame class adding some functionalities for easier Table manipulation.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @classmethod
    def from_csv(cls, filepath_or_buffer, **kwargs):
        """
        Read Table data from a csv file. This class method actually wraps around pandas.read_csv() function in order to
         return a Table object instead of a DataFrame. So any parameter that can be passed to pandas.read_csv() is also
         available for this Table.from_csv()

        :param filepath_or_buffer: string, path object or file-like object (see pandas.read_csv).
         Any valid path string is acceptable, including a valid URL (http, ftp, s3, gs, file).
         Path objects os.PathLike are acceptable.
         File-like objects are any object with a read() method, such as file handlers or StringIO.
        :param kwargs: any extra parameters to pass to pandas.read_csv() function.
        :return: a Table object.
        """
        return Table(data=read_csv(filepath_or_buffer, **kwargs))

    def set_index(self, keys, drop=True, append=False, inplace=False, verify_integrity=False):
        """
        Overwrite DataFrame.set_index method in order to return a Table object.

        :param keys: the index key(s) either as a single label or an array of labels (Series, Index, numpy.ndarray or an
         instance of Iterator
        :param drop: bool, default True
         Delete columns to be used as the new index.
        :param append: bool, default False
         Whether the specified columns should be append to existing index.
        :param inplace: bool, default False
         Modify the Table in place instead of creating a new instance.
        :param verify_integrity: bool, default False
         Check the new index for duplicates.
        :return: Table
         With a new index.
        """
        return Table(data=super().set_index(keys, drop, append, inplace, verify_integrity))

    def get_common_keys(self, other):
        """
        Return the index keys in common between self Table and other Table

        :param other: the other Table object.
        :return: Table
         The keys that are common between the two tables.
        """
        return Table(data=[x for x in self.index if x in other.index], columns=self.index.names)

    def get_keys_not_in(self, other):
        """
        Return the index keys in self Table which are not in other Table

        :param other: the other Table object.
        :return: Table
         The keys that are in self Table but not in the other one.
        """
        return Table(data=[x for x in self.index if x not in other.index], columns=self.index.names)
