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

    def join(self, other, on=None, how="left", lsuffix="", rsuffix="", sort=False):
        """
        Overwrite DataFrame.join method in order to return a Table object.

        :param other: the other Table (or DataFrame) to join
        :param on: str, list of str or array-like
         Column or index level name(s) on the caller to join on the index of other. If None, then joins index on index.
        :param how: {'left', 'right', 'outer', 'inner'}, default 'left'
         How to handle the operation of the two objects. See DataFrame.join() documentation for more details.
        :param lsuffix: str, default ''
         Suffix to use for left frames's overlapping columns
        :param rsuffix: str, default ''
         Suffix to use for right frames's overlapping columns
        :param sort: bool, default False
         Order result DataFrame lexicographically by the join key. If False, the order of the join key depends on the
        join type (how keyword).
        :return: A Table containing columns from both the caller and `other`
        """
        return Table(data=super().join(other, on, how, lsuffix, rsuffix, sort))

    def get_common_keys(self, other, index_self=None, index_other=None):
        """
        Return the index keys in common between self Table and other Table

        :param other: the other Table object
        :param index_self: a list of references to the columns defining indices in self Table. None refers
         to the original index
        :param index_other: a list of references to the columns defining indices in other Table. None refers
         to the original index
        :return: Table
         The keys that are common between the two tables.
        """
        ego = self if index_self is None else self.set_index(index_self)
        alter = other if index_other is None else other.set_index(index_other)
        return Table(data=[x for x in ego.index if x in alter.index], columns=ego.index.names)

    def get_common_rows(self, other, index_self=None, index_other=None, lsuffix="", rsuffix="", drop=True):
        """
        Return rows with the same index values. Index columns from the other Table object are droppped by default but
         can be kept if drop is set to False.

        :param other: the other Table object
        :param index_self: a list of references to the columns defining indices in self Table. None refers
         to the original index
        :param index_other: a list of references to the columns defining indices in other Table. None refers
         to the original index
        :param lsuffix: str,
         Suffix to add to left column names in case of redundancy.
        :param rsuffix: str,
         Suffix to add to right column names in case of redundancy.
        :param drop: True by default. Drop the other index columns in the resulting Table.
        :return: Table with the same index as self Table
         The rows that have the same key values between the two tables.
        """
        ego = self.copy() if index_self is None else self.set_index(index_self)
        alter = other.copy() if index_other is None else other.set_index(index_other, drop=drop)
        alter.index.names = ego.index.names
        return ego.join(alter, how='inner', lsuffix=lsuffix, rsuffix=rsuffix).reset_index()

    def get_keys_not_in(self, other, index_self=None, index_other=None):
        """
        Return the index keys in self Table which are not in other Table

        :param other: the other Table object
        :param index_self: a list of references to the columns defining indices in self Table. None refers
         to the original index
        :param index_other: a list of references to the columns defining indices in other Table. None refers
         to the original index
        :return: Table
         The keys that are in self Table but not in the other one.
        """
        ego = self.copy() if index_self is None else self.set_index(index_self, drop=False)
        alter = other.copy() if index_other is None else other.set_index(index_other, drop=False)
        return Table(data=[x for x in ego.index if x not in alter.index], columns=ego.index.names)
