"""
Extension module to the Biopython Bio.SeqUtils module
"""
import re

from Bio.Data import IUPACData
from gffpandas.gffpandas import Gff3DataFrame
from pandas import DataFrame


def ambiguous2string(code, protein=False):
    """
    Converts an ambiguous residue into a string with all compatible unambiguous residues. If the input code is not
    ambiguous, it is returned without any conversion.
    :param code: the input code to be converted into a list of residues.
    :param protein: True if residue is amino-acid
    :return: a string corresponding to the unambiguous residues compatible with the input code
    """
    if protein is True:
        if code.upper() in 'ACDEFGHIKLMNPQRSTVWWY':
            return code.upper()
        if code.upper() in 'BZX':
            return {'B': 'DN', 'Z': 'EQ', 'X': '.'}[code.upper()]
        raise ValueError('Invalid amino-acid %s' % code)
    if code.upper() in IUPACData.ambiguous_dna_values:
        return ''.join(sorted(IUPACData.ambiguous_dna_values[code.upper()]))
    raise ValueError('Invalid nucleotide %s' % code)


def isambiguous(code, protein=False):
    """
    Checks code is an ambiguous residue specification or not.
    :param code: the input code that must be checked for ambiguity
    :param protein: True if code is amino-acid code
    :return: Boolean True if code is ambiguous, False otherwise
    """
    return (protein is False and code in 'RYSWKMBDHVN') or (protein is True and code in 'BZX')


def pattern2regex(pattern, protein=False):
    """
    Converts a fuzznuc or fuzzpro-like pattern into a regular expression that can be used to search a sequence string.
    [ABC] => any of ABC residues
    {ABC} => any residue except ABC
    <ABC... => start of sequence
    ...ABC> => end of sequence
    A(n)(ABC)(n) => repeat residue or subsequence n times
    A(n,m)(ABC)(n,m) => repeat residue or subsequence from n up to m times

    :param pattern: the pattern definition (string)
    :param protein: True if pattern applies to a protein sequence, False otherwise.
    :return: the regular expression pattern as a string
    """
    pattern = re.sub(r'{(\w+)}', r'[^\1]', pattern)  # none of the specified residues
    pattern = re.sub(r'^<', r'^', pattern)  # start of sequence
    pattern = re.sub(r'>$', r'$', pattern)  # end of sequence
    pattern = re.sub(r'\((\d+(,\d*)*)\)', r'{\1}', pattern)  # repeat residue or subsequence

    # replace ambiguous code by the corresponding list of residues, adding [] if ot already within [], to produce
    # the regular expression
    regex = ''
    flag = False
    for c in pattern:
        if isambiguous(c, protein=protein) is True:
            c = ambiguous2string(c, protein=protein)
            if flag is False and c != '.':
                c = '[' + c + ']'
        elif c == '[':
            flag = True
        elif c == ']':
            flag = False
        regex += c

    return regex


class SeqFilter:
    """
    Class for the creation of a sequence filter to specify filtering criteria and applying the filter to a list of
    sequence records.
    minlength: minimum length of sequence
    maxlength: maximum length of sequence
    pattern: sequence pattern
    name: sequence name
    keep: boolean, if True, keep the records respecting the criteria, otherwise, discard them and keep the others
    """

    def __init__(self):
        self._minlength = None
        self._maxlength = None
        self._pattern = None
        self._name = None
        self._keep = True

    def length(self, minlength=None, maxlength=None):
        """
        Minimal and maximal length specification
        :param minlength: minimal accepted length
        :param maxlength: maximal accepted length
        :return: SeqFilter instance
        """
        self._minlength = minlength
        self._maxlength = maxlength
        return self

    def pattern(self, pattern=None):
        """
        pattern specification
        :param pattern: pattern that must be in the sequence
        :return: SeqFilter instance
        """
        self._pattern = pattern
        return self

    def name(self, name=None):
        """
        sequence record name specification
        :param name: name regular expression
        :return: SeqFilter instance
        """
        self._name = name
        return self

    def keep(self, keep=True):
        """
        Boolean defining whether the matching sequences must be kept (True) or removed (False)
        :param keep: True to keep positive sequences, False to remove them
        :return: SeqFilter instance
        """
        self._keep = keep
        return self

    def length_applies(self, r):
        """
        test whether length criterion applies to the sequence record
        :param r: the sequence record to test
        :return: boolean True if criterion applies or False otherwise
        """
        if self._keep is True:
            return r.seq.length_in_range(self._minlength, self._maxlength)
        return (self._minlength, self._maxlength) == (None, None) or r.seq.length_in_range(self._minlength,
                                                                                           self._maxlength) is False

    def pattern_applies(self, r):
        """
        test whether parameter criterion applies to the sequence record
        :param r: the sequence record to test
        :return: boolean True if criterion applies or False otherwise
        """
        if self._keep is True:
            return (self._pattern is None) or (r.seq.search(self._pattern) != [])
        return (self._pattern is None) or (r.seq.search(self._pattern) == [])

    def name_applies(self, r):
        """
        test whether name criterion applies to the sequence record
        :param r: the sequence record to test
        :return: boolean True if criterion applies or False otherwise
        """
        if self._keep is True:
            return (self._name is None) or (re.match(self._name, r.name) is not None)
        return (self._name is None) or (re.match(self._name, r.name) is None)

    def apply(self, records):
        """
        Filters a list of SeqRecords instances, keeping only records satisfying the specified criteria of length, match
        of a pattern, name specification. It is possible to invert the filtering process by setting the keep boolean to
        False and thus only keep records which do not satisfy the criteria.
        :param records: list of SeqRecord instances to apply
        """
        filtered = []
        for r in records:
            if all([self.length_applies(r),
                    self.pattern_applies(r),
                    self.name_applies(r)
                    ]):
                filtered.append(r)
        return filtered


class GFF(Gff3DataFrame):
    """
    Manipulation of features based upon gffpandas package
    """
    def __init__(self, feature_list=None):
        super().__init__(input_df=self.df_from_feature_list(feature_list))

    @classmethod
    def df_from_feature_list(cls, feature_list):
        """
        Create a pandas DataFrame from a list of features 'SeqFeatureEM2 or SeqFeature)
        :param feature_list: the list of features
        :return: the dataframe object
        """
        df = DataFrame(columns=['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
        if feature_list is not None:
            for ft in feature_list:
                df = df.append(cls.df_from_feature(ft), ignore_index=True)
        return df

    @staticmethod
    def df_from_feature(ft):
        """
        Create a pandas DataFrame from a feature (SeqFeatureEM2 or SeqFeature)
        :param ft: the feature to convert into a dataframe
        :return: the resulting dataframe
        """
        source = ft.qualifiers['source'] if 'source' in ft.qualifiers else ''
        score = ft.qualifiers['score'] if 'score' in ft.qualifiers else '0'
        if 'phase' in ft.qualifiers:
            phase = ft.qualifiers['phase']
        elif 'frame' in ft.qualifiers:
            phase = ft.qualifiers['frame']
        else:
            phase = '0'
        strand = ['-', '?', '+'][ft.location.strand + 1]
        return DataFrame([[ft.location.ref, source, ft.type, str(int(ft.location.start)), str(int(ft.location.end)),
                           score, strand, phase, ';'.join([k + '=' + v for k, v in ft.qualifiers.items()])]],
                         columns=['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
