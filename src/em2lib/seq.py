"""
Extension of Bio.Seq.Seq class from Biopython to add or improve functionalities
"""
import re
from typing import List, Union, Match

from Bio.Seq import Seq
from pandas import DataFrame
from pandas import concat

import em2lib.seq_utils


class SeqEM2(Seq):
    """
    SeqEM2 class providing extension to Bio.Seq.Seq class of BioPython package.
    """

    def __init__(self, data, seqtype):
        super().__init__(data)
        self.seqtype = seqtype
        self.any_residue = 'X' if seqtype == 'prot' else 'N'

    @classmethod
    def dna(cls, data):
        """
        Creates a DNA sequence

        :param data: the sequence string
        :return: a SeqEM2 DNA instance
        """
        return cls(data, 'dna')

    @classmethod
    def protein(cls, data):
        """
        Creates a protein sequence

        :param data: the sequence string
        :return: a SeqEM2 protein instance
        """
        return cls(data, 'prot')

    def is_protein(self):
        """
        Tests whether sequence was created as a protein

        :return: boolean, True if sequence was created as a protein, False otherwise.
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

        :param pattern: the pattern to be searched for that is converted into a regular expression
        :return: a list of re.Match objects
        """
        regex = em2lib.seq_utils.pattern2regex(pattern, protein=self.is_protein())
        return self.re_search(regex)

    def get_orfs(self, start=['ATG'], stop=None):
        """
        Determines all open reading frames in a sequence. It only examines the forward strand. All
        the returned ORFs have a length that is a multiple of 3. Thus, for sequences without any
        stop codon, 3 ORFs are returned, one for each frame.

        :param start: list of accepted start codon or None if ORFs do not need to start at a start
         codon.
        :param stop: list of accepted stop codons
        :return: a set of tuples (start, end) of orfs where start is the starting positon of the
         orf and end its ending position, not including the stop codon
        """
        seq_end = DataFrame([[self.__len__() - 2, self.__len__(), (self.__len__() - 2) % 3],
                             [self.__len__() - 1, self.__len__(), (self.__len__() - 1) % 3],
                             [self.__len__(), self.__len__(), self.__len__() % 3]],
                            columns=['start', 'end', 'frame'])
        if stop is None:
            stop = ['TAG', 'TGA', 'TAA']
        stop_matches = self.re_search('|'.join(stop))
        if stop_matches:
            stop_codons = concat(
                [concat([DataFrame([[match.start(), match.end(), match.start() % 3]],
                                   columns=['start', 'end', 'frame']) for match in stop_matches]
                        , ignore_index=True),
                 seq_end], ignore_index=True)
        else:
            stop_codons = seq_end

        if start is not None:
            start_matches = self.re_search('|'.join(start))
            if start_matches:
                start_codons = concat([DataFrame([[match.start(), match.end(), match.start() % 3]],
                                                 columns=['start', 'end', 'frame'])
                                       for match in start_matches], ignore_index=True)
            else:
                start_codons = DataFrame([], columns=['start', 'end', 'frame'])
        else:
            seq_start = DataFrame([start_pos for start_pos in [[0, 3, 0], [1, 4, 1], [2, 5, 2]]
                                   if start_pos not in [[row.start, row.end, row.frame] for row in
                                                        stop_codons.itertuples()]],
                                  columns=['start', 'end', 'frame'])
            start_codons = concat([seq_start, concat([DataFrame([[row.end, row.end + 3, row.frame]],
                                                                columns=['start', 'end', 'frame'])
                                                      for row in stop_codons.itertuples()],
                                                     ignore_index=True)], ignore_index=True)
        orfs = set()
        for start_codon in start_codons.itertuples():
            for stop_codon in stop_codons.itertuples():
                if start_codon.frame == stop_codon.frame and stop_codon.start > start_codon.end:
                    orfs.add(tuple([start_codon.start, stop_codon.start]))
                    break
        return orfs
