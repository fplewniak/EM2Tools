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
