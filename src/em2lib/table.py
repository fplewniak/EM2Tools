"""
Some utilities for Table manipulation, comparison, etc. using pandas.DataFrame module.
"""
#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
from pandas import DataFrame


class Table:
    """
    A class adding some utilities for easier DataFrame manipulation, comparison, etc.
    """

    @staticmethod
    def get_common_keys(first, other, keys_first=None, keys_other=None):
        """
        Return the index keys in common between first DataFrame and other DataFrame

        :param first: the first DataFrame object to compare
        :param other: the other DataFrame object to compare
        :param keys_first: a list of references to the columns defining indices in first DataFrame. None refers
         to the original index
        :param keys_other: a list of references to the columns defining indices in other DataFrame. None refers
         to the original index
        :return: DataFrame, the keys that are common between the two DataFrames.
        """
        ego = first if keys_first is None else first.set_index(keys_first)
        alter = other if keys_other is None else other.set_index(keys_other)
        return DataFrame(data=[x for x in ego.index if x in alter.index], columns=ego.index.names)

    @staticmethod
    def get_common_rows(first, other, keys_first=None, keys_other=None, lsuffix="", rsuffix="", drop=True):
        """
        Return rows with the same key values. Key columns from the other DataFrame object are droppped by default
         but can be kept if drop is set to False.

        :param first: the first DataFrame object to compare
        :param other: the other DataFrame object to compare
        :param keys_first: a list of references to the columns defining keys in first DataFrame. None refers
         to the original index
        :param keys_other: a list of references to the columns defining keys in other DataFrame. None refers
         to the original index
        :param lsuffix: str,
         Suffix to add to first column names in case of redundancy.
        :param rsuffix: str,
         Suffix to add to right column names in case of redundancy.
        :param drop: True by default. Drop the other key columns in the resulting DataFrame.
        :return: DataFrame with the same keys as first DataFrame, containing the rows that have the same key values
         between the two DataFrames.
        """
        ego = first.copy() if keys_first is None else first.set_index(keys_first)
        alter = other.copy() if keys_other is None else other.set_index(keys_other, drop=drop)
        alter.index.names = ego.index.names
        return ego.join(alter, how='inner', lsuffix=lsuffix, rsuffix=rsuffix).reset_index()

    @staticmethod
    def get_keys_not_in(first, other, keys_first=None, keys_other=None):
        """
        Return the keys in first DataFrame which are not in other DataFrame

        :param first: the first DataFrame object to compare
        :param other: the other DataFrame object to compare
        :param keys_first: a list of references to the columns defining keys in first DataFrame. None refers
         to the original index
        :param keys_other: a list of references to the columns defining keys in other DataFrame. None refers
         to the original index
        :return: DataFrame, the keys that are in first DataFrame but not in the other one.
        """
        ego = first.copy() if keys_first is None else first.set_index(keys_first, drop=False)
        alter = other.copy() if keys_other is None else other.set_index(keys_other, drop=False)
        return DataFrame(data=[x for x in ego.index if x not in alter.index], columns=ego.index.names)

    @staticmethod
    def get_rows_not_in(first, other, keys_first=None, keys_other=None):
        """
        Return a DataFrame with rows with keys that are in first but not in other.

        :param first: the first DataFrame object to compare
        :param other: the other DataFrame object to compare
        :param keys_first: a list of references to the columns defining keys in first DataFrame. None refers
         to the original index
        :param keys_other: a list of references to the columns defining keys in other DataFrame. None refers
         to the original index
        :return: DataFrame, the rows whose keys are in first DataFrame but not in the other one.
        """
        ego = first.copy() if keys_first is None else first.set_index(keys_first)
        alter = other.copy() if keys_other is None else other.set_index(keys_other)
        return DataFrame(data=ego.loc[[x for x in ego.index if x not in alter.index]]).reset_index()
