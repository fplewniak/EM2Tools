"""
Some utilities for Table manipulation, comparison, etc. using pandas.DataFrame module.
"""
#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
from pandas import DataFrame
from pandas import Series
import numpy


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


class TableTransform():
    """
    Class of objects to transform a DataFrame according to different criteria. Different independent transformations
     can be chained since all transformation methods return self object.
    """

    def __init__(self, df):
        """
        TableTransform object creator

        :param df: The original DataFrame to transform
        """
        self.orgnl_df = df
        self.wrkg_df = df.copy()

    def result(self):
        """
        Getter method to retrieve the DataFrame resulting from the transformations
        :return: the transformed DataFrame
        """
        return self.wrkg_df

    def update(self):
        """
        Update orginal DataFrame with the current working copy. This enables the application of conditions to previously
         modified values. This action cannot be undone, so it is recommended to use it after all modifications based on
         original data have been performed.
        """
        self.orgnl_df = self.wrkg_df
        return self

    def cond_transform(self, cond=lambda x: True, iftrue=lambda x: x, columns=None, original=True):
        """
        Conditional transform of a DataFrame. If condition function returns True, then the iftrue function is applied
         to transform the corresponding element.

        :param cond: function, DataFrame or Series
         A function is applied to elements of the original DataFrame in order to define which elements in the
          working DataFrame should be transformed.
         A DataFrame is a mask of bool values with the same shape as the original DataFrame. Elements in the working
         A Series is a mask of bool values defining in which rows the transformation should be applied.
        :param iftrue: function or a DataFrame with the same shape as te original table.
         the function to apply if condition is True, returns a single value from a single value. This
         function should test whether it is applicable to its input and return the original value if not.
        :param columns: str or list thereof specifying the column(s) which the transformation should be applied to
        :param original: if True, the original DataFrame is used to assess the condition function, otherwise, the
         working copy is used. This has an effect only if cond is a function and therefore only works for elementwise
         transformation. For rowwise selection, cond should be a DataFrame that can be computed by a function taking as
         input the current result DataFrame self.result()
        :return: this TableTransform instance
        """
        if isinstance(cond, (DataFrame, Series)):
            mask = cond
        elif original:
            mask = self.orgnl_df.applymap(cond)
        else:
            mask = self.wrkg_df.applymap(cond)
        if isinstance(iftrue, DataFrame):
            tmp_df = self.wrkg_df.mask(mask, iftrue)
        else:
            tmp_df = self.wrkg_df.mask(mask, self.wrkg_df.applymap(iftrue))
        if columns is not None:
            self.wrkg_df.loc[:, columns] = tmp_df.loc[:, columns]
        else:
            self.wrkg_df = tmp_df
        return self

    def replace(self, columns=None, **kwargs):
        """
        A method wrapping the DataFrame.replace() method and adding the columns parameter to apply the replacement only
         in the specified columns.

        :param columns: str or list thereof specifying the column(s) to apply replacement to
        :param kwargs: named arguments to pass to DataFrame.replace() method. Note that inplace argument is deactivated
         as it is not needed here and might actually interfere with the columns argument. As a matter of fact, the
         working DataFrame is always modified.
        :return: this TableTransform instance
        """
        kwargs['inplace'] = False
        tmp_df = self.wrkg_df.replace(**kwargs)
        if columns is not None:
            self.wrkg_df.loc[:, columns] = tmp_df.loc[:, columns]
        else:
            self.wrkg_df = tmp_df
        return self

    def randomize(self, by=None, replacement=False, seed=None):
        """
        Randomize a DataFrame by row, column or both with or without replacement. If by is row (or column), then whole
         rows (or whole columns) are randomized. Otherwise, all elements of the DataFrame are resampled.

        :param by: specifies whether randomization should be performed by row, by column or by element (None)
        :param replacement: if True, then randomization will occur with replacement.
        :return: the current TableTransform object
        """
        tmp_df = DataFrame(columns=self.wrkg_df.columns)
        if by == 'column':
            for column in self.wrkg_df.columns:
                tmp_df[column] = self.wrkg_df[column].sample(frac=1, axis=0, replace=replacement,
                                                             random_state=seed).reset_index(drop=True)
        elif by == 'row':
            pass
        else:
            pass
        return self
