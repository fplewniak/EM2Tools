"""
Basic utilities that can be shared between different modules of the API
"""
#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
from numbers import Number


def norm_sum_to_one(data):
    """
    Normalization of values in a 1-D array so they sum to one

    :param data: the 1-D array to normalize
    :return: the normalized array
    """
    return data / data.sum() if all(isinstance(x, Number) for x in data) else data


def norm_max(data):
    """
    Normalization of values in a 1-D array relative to the max value

    :param data: the 1-D array to normalize
    :return: the normalized array
    """
    return data / data.max() if all(isinstance(x, Number) for x in data) else data


def norm_min_max(data):
    """
    Normalization of values in a 1-D array to make them fall between 0 (min value) and 1 (max value)

    :param data: the 1-D array to normalize
    :return: the normalized array
    """
    return (data - data.min()) / (data.max() - data.min()) if all(isinstance(x, Number) for x in data) else data


def zscore(data):
    """
    Computes the Z-score of values in a 1-D array

    :param data: the 1-D array to normalize
    :return: the normalized array
    """
    return (data - data.mean()) / data.std() if all(isinstance(x, Number) for x in data) else data
