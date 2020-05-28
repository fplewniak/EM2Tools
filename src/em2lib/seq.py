"""
Extension of Bio.Seq.Seq class from Biopython to add or improve functionalities
"""
import re
from typing import List, Union, Match

from Bio.Seq import Seq

import em2lib.seq_utils


class SeqEM2(Seq):
    """
    SeqEM2 class providing extension to Bio.Seq.Seq class of BioPython package.
    """

    def __init__(self, data, seqtype):
        super().__init__(data)
        self.seqtype = seqtype

    @classmethod
    def dna(cls, data):
        """
        Creates a DNA sequence

        :param data: The sequence string
        :return: a SeqEM2 DNA instance
        """
        return cls(data, 'dna')

    @classmethod
    def protein(cls, data):
        """
        Creates a protein sequence

        :param data: The sequence string
        :return: a SeqEM2 protein instance
        """
        return cls(data, 'prot')

    def is_protein(self):
        """
        Tests whether sequence was created as a protein

        :return:
        """
        return self.seqtype == 'prot'

    def length_in_range(self, minlength=None, maxlength=None):
        """
        Checks whether the sequence length is with the specified range.

        :param minlength: lower length bound
        :param maxlength: upper length bound
        :return: True if sequence length is within specified range, False otherwise
        """
        if minlength is None and maxlength is None:
            return True
        if minlength is None:
            return self.__len__() <= maxlength
        if maxlength is None:
            return minlength <= self.__len__()
        return minlength <= self.__len__() <= maxlength

    def re_search(self, regex):
        """
        Searches a sequence using a regular expression

        :param regex: the regular expression
        :return: a list of re.Match instances
        """
        result: List[Match[Union[str, bytes]]] = []
        for match in re.finditer(regex, self._data):
            result.append(match)
        return result

    def search(self, pattern):
        """
        Searches sequence for a pattern specified with a fuzznuc or fuzzpro like syntax

        :param pattern: the pattern to be searched for that is converted into a rgular expression
        :return: a list of re.Match objects
        """
        regex = em2lib.seq_utils.pattern2regex(pattern, protein=self.is_protein())
        return self.re_search(regex)
