"""
Some utilities for Table manipulation, comparison, etc. using pandas.DataFrame module.
"""
#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
from numbers import Number

from pandas import DataFrame
from pandas import Series
from pandas import MultiIndex
import numpy
import em2lib.utils as em2


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

    @staticmethod
    def statistics(table, groupby=None, columns=None, func=numpy.mean):
        """
        Returns a DataFrame containing the requested statistics on the specified columns in a table, optionally grouping
        data according to one or several other columns.

        :param table: the table to compute the statistics
        :param groupby: the column or list of columns for grouping data
        :param columns: the column or list of columns to compute the statistics about
        :param func: the statistical function or list thereof to be computed
        :return: a DataFrame with the requested statistics
        """
        if not isinstance(func, list):
            func = [func]
        if columns is None:
            columns = table.columns
        if not isinstance(columns, list):
            columns = [columns]
        if groupby is None:
            return table[columns].agg(func)
        stat_index = MultiIndex.from_product([[f.__name__ for f in func], columns])
        stat_df = DataFrame(columns=stat_index)
        for each_func in func:
            stat_df[[x for x in stat_index if x[0] == each_func.__name__]] = \
                table.groupby(groupby)[columns].apply(func=each_func)
        return stat_df

    @staticmethod
    def implode(input_df, index=None):
        """
        Reverts action of DataFrame.explode method. This method can also be applied to DataFrames that are not the
        result of an explode() call but in this case, applying back the explode() method may not yield the original
        DataFrame.

        :param input_df: the input DataFrame
        :param index: the index column(s)
        :return: the resulting DataFrame
        """
        # copy input DataFrame to avoid any modification of the original object
        exp_df = input_df.copy()
        # set the index if needed
        if index is not None:
            exp_df.set_index(index, drop=False, inplace=True)
        # prepare the DataFrame for output
        tmp_df = DataFrame(index=numpy.unique(exp_df.index), columns=exp_df.columns)
        # gather elements of each column in a list for each unique value in index
        for col in exp_df.columns:
            unique_lists = {}
            for i in exp_df.index:
                unique_lists[i] = []
                # if cell contains only one tuple element, then place it in a list to keep it as a single element
                cell = [exp_df.loc[i, col]] if isinstance(exp_df.loc[i, col], tuple) else exp_df.loc[i, col]
                # get unique elements in cell for the current index level
                for cell_element in cell:
                    if cell_element not in unique_lists[i]:
                        unique_lists[i].append(cell_element)
            # update the current column with the lists of unique elements in current cell
            tmp_df[col] = Series(unique_lists)
        # replace list with only one element by the element itself
        tmp_df = TableTransform(tmp_df).cond_transform(cond=lambda x: len(x) == 1, iftrue=lambda x: x[0]).result()
        # if a column index has been specified then reset the index to avoid duplication of data...
        if index is not None:
            return tmp_df.reset_index(drop=True)
        # ...otherwise, keep the index
        return tmp_df

    @staticmethod
    def collapse(input_df, groupby=None, columns=None, name=None):
        """
        Collapses a DataFrame into a single column containing lists of tuples representing all the values in the
        specified columns with the same index as defined by groupby.
        Note: results may be inconsistent if input DataFrame contains elements that are lists or tuples.

        :param input_df: the input DataFrame
        :param groupby: the column(s) defining the index for collapsing values
        :param columns: the columns that will be collapsed
        :param name: the name for the resulting column
        :return: a single-column collapsed DataFrame
        """
        # copy input DataFrame to avoid any modification of the original object
        exp_df = input_df.copy()
        # set the index if needed
        if groupby is not None:
            exp_df.set_index(groupby, drop=True, inplace=True)
        # if no column is specified, then use all of them except the index
        if columns is None:
            columns = list(exp_df.columns) if len(exp_df.columns) > 1 else exp_df.columns[0]
        # prepare the DataFrame for output
        if name is None:
            name = tuple(columns) if isinstance(columns, list) else columns
        tmp_df = DataFrame(columns=[name])
        # gather all elements of rows in lists of tuples, a tuple, or as a single element
        row_lists = {}
        for i in exp_df.index:
            # there are several values in multiple columns for the same index value, get a list of tuples
            if isinstance(exp_df.loc[i, columns], DataFrame):
                row_lists[i] = [tuple(row[1]) for row in exp_df.loc[i, columns].iterrows()]
            # there are several values in one column (list of values) or one value in several columns (tuple of values)
            elif isinstance(exp_df.loc[i, columns], Series):
                row_lists[i] = [row[1] for row in exp_df.loc[i, columns].items()]
                # there is one value in several columns, get a tuple of values
                if isinstance(columns, list):
                    row_lists[i] = tuple(row_lists[i])
            # there is one value in one column, get a single value
            else:
                row_lists[i] = [exp_df.loc[i, columns]] if isinstance(exp_df.loc[i, columns], list) \
                    else exp_df.loc[i, columns]
        tmp_df[name] = Series(row_lists)
        # replace list with only one element by the element itself and return the resulting DataFrame with reset index
        tmp_df[name] = TableTransform(DataFrame(tmp_df[name])).cond_transform(cond=lambda x: len(x) == 1,
                                                                              iftrue=lambda x: x[0]).result()
        tmp_df.index.rename(groupby, inplace=True)
        return tmp_df

    @staticmethod
    def expand(input_df, columns=None):
        """
        Expands a collapsed table. Exact reverse of Table.collapse() method.

        :param input_df: the collapsed table to expand
        :param columns: names of the columns, if None, will attempt to deduce them from the collapsed column name
         expected to be a tuple of the original names as produced by Table.collapse(). Note that if less column names
         are specified than the actual number of columns to expand then extra column names will be added automatically.
        :return: the expanded table
        """

        rows = []
        # get the list of elements for each value of the index
        for i in input_df.index:
            element_list = input_df.loc[i].to_list()
            # if the list of elements is a list of list of tuples, then keep only the list of tuples
            if len(numpy.array(element_list).shape) > 2:
                element_list = element_list[0]
            # append each list of elements and index values to the list of rows
            for element in element_list:
                rows.append(list(i) + list(element))
        # if no column names are specified then try to determine a list of names from the collapsed column name
        if columns is None:
            columns = list(input_df.columns)[0] if isinstance(input_df.columns[0], tuple) else list(input_df.columns)
        # if only a single name was passed then place it into a list
        elif not isinstance(columns, list):
            columns = [columns]
        # if there are less specified columns than data columns then add extra columns
        if len(rows[0]) > len(input_df.index.names + columns):
            columns = columns + ['col_' + str(x) for x in range(len(rows[0]) - len(input_df.index.names + columns))]
        # return the DataFrame constructed from the list of rows
        return DataFrame(rows, columns=input_df.index.names + columns)

    @staticmethod
    def table_from_edges(df_edges: DataFrame, graph=True, directed=True, link=None, no_link=None):
        """
        Turns a DataFrame representing edges in two or three columns into a table representing edges between nodes in
         rows and columns. The first two columns of the input DataFrame represent the linked nodes. The optional third
         column contains a value (edge weight, label, etc.).
         A call to table_from_edges(df_edges, no_link=0, link=1, graph=False) is basically equivalent to a call to
         pandas.crosstab(df_edges[0], df_edges[1]).
        A call to table_from_edges(df_edges, graph=False) is basically equivalent to a call to df_edges.pivot(0,1,2).

        :param df_edges: the edge DataFrame,
        :param graph: if True, the edge list represents a graph and the resulting table will be a square: all nodes
         in the first two columns of the edge DataFrame are placed in both the index and the columns of the graph table.
         If False, the edge list represents a relation between two sets. The first column of the edge DataFrame will
         define rows (first set) and the second column will define the columns (second set).
        :param directed: specifies if the graph is directed or not. If True, the resulting table will be symmetrical.
         This parameter has no effect if graph is False
        :param link: the value to associate the specified edges with. This value will replace the original value in
         df_edges DataFrame or specify one for two-column edge DataFrames. This value is needed for two-column edge
         definition.
        :param no_link: the value to associate the unspecified edges with. NaN will be used if this parameter is not
         specified.
        :return: a table DataFrame
        """
        if graph:
            # create a square table with all nodes in df_edges
            nodes = (df_edges[0].append(df_edges[1])).unique()
            df_table = DataFrame(index=nodes, columns=nodes)
        else:
            # create an asymmetrical table with row nodes and column nodes in df_edges
            df_table = DataFrame(index=df_edges[0].unique(), columns=df_edges[1].unique())
        # populate the table
        for edge in df_edges.iterrows():
            if link is not None:
                # set the specified edge with the value specified in link parameter
                df_table.loc[edge[1].loc[0], edge[1].loc[1]] = link
                if graph is True and directed is False:
                    # in the case of a directed graph, set also the symmetrical edge
                    df_table.loc[edge[1].loc[1], edge[1].loc[0]] = link
            else:
                # set the specified edge with the value specified in the third column of df_edges
                df_table.loc[edge[1].loc[0], edge[1].loc[1]] = edge[1].loc[2]
                if graph is True and directed is False:
                    # in the case of a directed graph, set also the symmetrical edge
                    df_table.loc[edge[1].loc[1], edge[1].loc[0]] = edge[1].loc[2]
        if no_link is not None:
            # replace NaN (unspecified edges) with the value specified in no_link parameter
            df_table.fillna(no_link, inplace=True)
        return df_table


class TableTransform:
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
        Getter method to retrieve the DataFrame resulting from the transformation

        :return: the transformed DataFrame
        """
        return self.wrkg_df

    def update(self):
        """
        Update original DataFrame with the current working copy. This enables the application of conditions to
        previously modified values. This action cannot be undone, so it is recommended to use it after all
        modifications based on original data have been performed.

        """
        self.orgnl_df = self.wrkg_df
        return self

    def cond_transform(self, cond=lambda x: True, iftrue=lambda x: x, columns=None, original=True):
        """
        Conditional transform of a DataFrame. If condition function returns True, then the iftrue function is applied
        to transform the corresponding element.

        :param cond: string, function, DataFrame or Series
         A string considered as a test function on the columns of the DataFrame working copy. This string is passed to
         self.wrkg_df.eval() for evaluation to produce a mask DataFrame.
         A function is applied to elements of the original DataFrame in order to define which elements in the
         working DataFrame should be transformed.
         A DataFrame is a mask of bool values with the same shape as the original DataFrame. Elements in the working
         A Series is a mask of bool values defining in which rows the transformation should be applied.
        :param iftrue: string, function or a DataFrame with the same shape as te original table.
         The function to apply if condition is True, returns a single value from a single value. This
         function should test whether it is applicable to its input and return the original value if not.
         If a DataFrame, it defines the values that will replace those where the condition is True.
         If a string, it is passed to self.wrkg_df.eval() for evaluation and is used as a DataFrame.
        :param columns: str or list thereof specifying the column(s) which the transformation should be applied to
        :param original: if True, the original DataFrame is used to assess the condition function, otherwise, the
         working copy is used. This has an effect only if cond is a function and therefore only works for elementwise
         transformation. For rowwise selection, cond should be a DataFrame that can be computed by a function taking as
         input the current result DataFrame self.result()
        :return: this TableTransform instance
        """
        if isinstance(cond, str):
            mask = self.wrkg_df.eval(cond)
        elif isinstance(cond, (DataFrame, Series)):
            mask = cond
        elif original:
            mask = self.orgnl_df.applymap(cond)
        else:
            mask = self.wrkg_df.applymap(cond)

        if isinstance(iftrue, str):
            iftrue = self.wrkg_df.eval(iftrue)

        if isinstance(iftrue, DataFrame):
            tmp_df = self.wrkg_df.mask(mask, iftrue)
        else:
            tmp_df = self.wrkg_df.mask(mask, self.wrkg_df.applymap(iftrue))

        if columns is not None:
            self.wrkg_df.loc[:, columns] = tmp_df.loc[:, columns]
        else:
            self.wrkg_df = tmp_df
        return self

    def combine(self, other, func, columns=None, **kwargs):
        """
        A method wrapping the DataFrame.combine() method and adding the columns parameter to apply the method only
        to the specified columns.

        :param other: the other DataFrame to combine to the current working DataFrame
        :param func: function that takes two series as inputs and returns a Series or a scalar to merge the two
         dataframes column by column.
        :param columns: str or list thereof specifying the column(s) which should be kept in the final result.
        :param kwargs: named arguments to pass to DataFrame.replace() method.
        :return: this TableTransform instance
        """
        if columns is not None:
            self.wrkg_df.loc[:, columns] = self.wrkg_df.loc[:, columns].combine(other.loc[:, columns], func, **kwargs)
        else:
            self.wrkg_df = self.wrkg_df.combine(other, func, **kwargs)
        return self

    def normalize(self, columns=None, norm=em2.norm_sum_to_one, by=None):
        """
        Normalize values in DataFrame, either by row, column or over all the values.

        :param columns: the columns whose values should be normalized
        :param norm: the normalization function, should apply to a 1-D array, default is ''norm_sum_to_one''
        :param by: defines whether the normalization should be done by column or by row. If None, then all values are
         used.
        :return: this TableTransform instance
        """
        # if no column is specified then use all of them...
        if columns is None:
            columns = self.wrkg_df.columns
        # ... but remove non-numerical columns
        for col in columns:
            if not all([isinstance(x, Number) for x in self.wrkg_df.loc[:, col]]):
                columns = columns.drop(col)
        # make a copy of the specified columns as a DataFrame for the normalization
        tmp_df = self.wrkg_df.loc[:, columns].copy()
        # if normalizing by row or by column,
        if by is not None:
            # translate to the appropriate axis parameter and apply the normalization on the requested axis
            axis = 0 if by == 'col' else 1
            tmp_df = DataFrame(tmp_df).apply(norm, axis=axis)
        else:
            # else, apply the normalization on all values
            tmp_df = DataFrame(numpy.reshape(DataFrame(tmp_df.to_numpy().flatten()).apply(norm).to_numpy(),
                                             tmp_df.shape), columns=tmp_df.columns)
        # send the normalized values in the original DataFrame
        self.wrkg_df.loc[:, tmp_df.columns] = tmp_df
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

    def randomize(self, by=None, replacement=False, columns=None, seed=None):
        """
        Randomize a DataFrame by row, column or all elements across the whole table with or without replacement.
        If by=row (or column), then elements of each row (or each column) are randomized within the row (or the column)
        Otherwise, all elements of the DataFrame are resampled across the whole table independently of rows or columns.

        :param by: specifies whether randomization should be performed by row, by column or by element (None)
        :param replacement: if True, then randomization will occur with replacement.
        :param columns: str or list thereof specifying the column(s) which the columns affected by the randomization
        :param seed: seed for random numbers generation passed to sample() method as the random_state argument
        :return: the current TableTransform object
        """
        if by == 'column':
            tmp_df = DataFrame(columns=self.wrkg_df.columns)
            for column in self.wrkg_df.columns:
                tmp_df[column] = self.wrkg_df[column].sample(frac=1, axis=0, replace=replacement,
                                                             random_state=seed).reset_index(drop=True)
        elif by == 'row':
            tmp_df = DataFrame(columns=self.wrkg_df.columns)
            for row in self.wrkg_df.index:
                rand_row = DataFrame(list(self.wrkg_df.loc[row, :])).transpose().sample(frac=1, axis=1,
                                                                                        replace=replacement,
                                                                                        random_state=seed)
                rand_row.columns = self.wrkg_df.columns
                tmp_df = tmp_df.append(rand_row)
            tmp_df = tmp_df.reset_index(drop=True)
        else:
            df_array = self.wrkg_df.to_numpy()
            tmp_df = DataFrame(numpy.reshape(
                numpy.random.choice(df_array.flat, size=len(df_array.flat), replace=replacement),
                df_array.shape), columns=self.wrkg_df.columns)

        if columns is not None:
            self.wrkg_df.loc[:, columns] = tmp_df.loc[:, columns]
        else:
            self.wrkg_df = tmp_df
        return self
