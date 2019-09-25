"""
Extension of Bio.Seq.Seq class from Biopython to add or improve functionalities
"""
import re
from typing import List, Union, Match

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import EM2Libs.SeqUtils


class SeqEM2(Seq):
    """
    SeqEM2 class providing extension to Bio.Seq.Seq class of BioPython package.
    """

    def __init__(self, data, alphabet):
        super().__init__(data, alphabet)

    @classmethod
    def dna(cls, data):
        """
        Creates a DNA sequence
        :param data: The sequence string
        :return: a SeqEM2 DNA instance
        """
        return cls(data, IUPAC.ExtendedIUPACDNA)

    @classmethod
    def protein(cls, data):
        """
        Creates a protein sequence
        :param data: The sequence string
        :return: a SeqEM2 protein instance
        """
        return cls(data, IUPAC.ExtendedIUPACProtein)

    def is_protein(self):
        """
        Tests whether sequence was created as a protein
        :return:
        """
        return self.alphabet == IUPAC.ExtendedIUPACProtein

    def re_search(self, regex):
        """
        Searches a sequence using a regular expression
        :param regex: the regular expression
        :return: a list of re.Match instances
        """
        result: List[Match[Union[str, bytes]]] = []
        for m in re.finditer(regex, self._data):
            result.append(m)
        return result

    def search(self, pattern):
        """
        Searches sequence for a pattern specified with a fuzznuc or fuzzpro like syntax
        :param pattern: the pattern to be searched for that is converted into a rgular expression
        :return: a list of re.Match objects
        """
        regex = EM2Libs.SeqUtils.pattern2regex(pattern, protein=self.is_protein())
        return self.re_search(regex)
