"""
Extension module to the Biopython Bio.SeqUtils module
"""
import re

from Bio.Data import IUPACData
from Bio.SeqFeature import FeatureLocation
from pandas import DataFrame
from gffpandas.gffpandas import Gff3DataFrame
from em2lib.seq_feature import SeqFeatureEM2


def ambiguous2string(code, protein=False):
    """
    Converts an ambiguous residue into a string with all compatible unambiguous residues.
    If the input code is not ambiguous, it is returned without any conversion.

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

    :return: boolean, True if code is ambiguous, False otherwise
    """
    return (protein is False and code in 'RYSWKMBDHVN') or (protein is True and code in 'BZX')


def pattern2regex(pattern, protein=False):
    """
    Converts a fuzznuc or fuzzpro-like pattern into a regular expression that can be used to search
    a sequence string.

    - [ABC] => any of ABC residues,
    - {ABC} => any residue except ABC,
    - <ABC... => start of sequence,
    - ...ABC> => end of sequence,
    - A(n)(ABC)(n) => repeat residue or subsequence n times,
    - A(n,m)(ABC)(n,m) => repeat residue or subsequence from n up to m times.

    :param pattern: the pattern definition (string)
    :param protein: True if pattern applies to a protein sequence, False otherwise.

    :return: the regular expression pattern as a string
    """
    pattern = re.sub(r'{(\w+)}', r'[^\1]', pattern)  # none of the specified residues
    pattern = re.sub(r'^<', r'^', pattern)  # start of sequence
    pattern = re.sub(r'>$', r'$', pattern)  # end of sequence
    pattern = re.sub(r'\((\d+(,\d*)*)\)', r'{\1}', pattern)  # repeat residue or subsequence

    # replace ambiguous code by the corresponding list of residues, adding [] if ot already
    # within [], to produce the regular expression
    regex = ''
    flag = False
    for code in pattern:
        if isambiguous(code, protein=protein) is True:
            code = ambiguous2string(code, protein=protein)
            if flag is False and code != '.':
                code = '[' + code + ']'
        elif code == '[':
            flag = True
        elif code == ']':
            flag = False
        regex += code

    return regex


class SeqFilter:
    """
    A class for the creation of a sequence filter to specify filtering criteria and applying the\
     filter to a list of sequence records.
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

    def length_applies(self, rec):
        """
        test whether length criterion applies to the sequence record

        :param rec: the sequence record to test
        :return: boolean True if criterion applies or False otherwise
        """
        if self._keep is True:
            return rec.seq.length_in_range(self._minlength, self._maxlength)
        return (self._minlength, self._maxlength) == (None, None) or rec.seq.length_in_range(
            self._minlength,
            self._maxlength) is False

    def pattern_applies(self, rec):
        """
        test whether parameter criterion applies to the sequence record

        :param rec: the sequence record to test
        :return: boolean True if criterion applies or False otherwise
        """
        if self._keep is True:
            return (self._pattern is None) or (rec.seq.search(self._pattern) != [])
        return (self._pattern is None) or (rec.seq.search(self._pattern) == [])

    def name_applies(self, rec):
        """
        test whether name criterion applies to the sequence record

        :param rec: the sequence record to test
        :return: boolean True if criterion applies or False otherwise
        """
        if self._keep is True:
            return (self._name is None) or (re.match(self._name, rec.name) is not None)
        return (self._name is None) or (re.match(self._name, rec.name) is None)

    def apply(self, records):
        """
        Filters a list of SeqRecords instances, keeping only records satisfying the specified
        criteria of length, match of a pattern, name specification. It is possible to invert
        the filtering process by setting the keep boolean to False and thus only keep records
        which do not satisfy the criteria.

        :param records: list of SeqRecord instances to apply
        :return: the filtered list of records
        """
        filtered = []
        for rec in records:
            if all([self.length_applies(rec),
                    self.pattern_applies(rec),
                    self.name_applies(rec)
                    ]):
                filtered.append(rec)
        return filtered


class GFF(Gff3DataFrame):
    """
    Manipulation of features based upon gffpandas package
    """

    def __init__(self, feature_list=None, input_df=None):
        init_df = DataFrame(
            columns=['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                     'attributes'])
        super().__init__(input_df=init_df)
        self.add_feature_list(feature_list)
        if input_df is not None:
            self.df = self.df.append(input_df)

    def add_feature_list(self, feature_list=None):
        """
        Adds a list of feature to the list of an existing GFF object

        :param feature_list: list of features to add to DataFrame
        :return: the GFF object with feature list appended
        """
        if feature_list is not None:
            for feature in feature_list:
                self.df = self.df.append(self.df_from_feature(feature), ignore_index=True)
        return self

    @staticmethod
    def df_from_feature(feature):
        """
        Create a pandas DataFrame from a feature (SeqFeatureEM2 or SeqFeature)

        :param feature: the feature to convert into a dataframe
        :return: the resulting dataframe
        """
        if feature is None:
            return DataFrame(
                columns=['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                         'attributes'])

        source = feature.qualifiers.get('source', 'unknown')
        score = feature.qualifiers.get('score', '0')
        phase = feature.qualifiers.get('phase', feature.qualifiers.get('frame', '0'))
        strand = ['-', '?', '+'][feature.location.strand + 1]
        qual_off = ['phase', 'source', 'score']
        qualifiers = ';'.join(
            [k + '=' + v for k, v in feature.qualifiers.items() if k not in qual_off]
            + ['id=' + feature.id])
        return DataFrame(
            [[feature.parent.id, source, feature.type, str(int(feature.location.start)),
              str(int(feature.location.end)),
              score, strand, phase, qualifiers]],
            columns=['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand',
                     'phase', 'attributes'])

    def to_feature_list(self, parents=None):
        """
        Converts features in a GFF object into a list of SeqFeatureEM2 objects

        :param parents: list of references to parent SeqRecord objects or a single parent reference
         if all features are defined in the same parent. If it is a list, it should be of the same
         length as the dataframe, repeating references as needed to get the right number.
        :return: a list of SeqFeatureEM2 objects
        """
        if isinstance(parents, list) is False:
            parents = [parents] * len(self.df)
        if len(parents) != len(self.df):
            raise ValueError(
                'The number of parents should match the number of features unless all features'
                ' have the same parent, in which case only one reference can be specified.')
        feature_list = []
        strands = {'-': -1, '?': 0, '+': 1}
        for index, row in self.df.iterrows():
            qualifiers = {q.split('=')[0]: q.split('=')[1] for q in row['attributes'].split(';')}
            for column in ['source', 'phase', 'score']:
                qualifiers[column] = row[column]
            feature_list.append(
                SeqFeatureEM2(parent=parents[index],
                              location=FeatureLocation(int(row['start']), int(row['end'])),
                              type=row['type'], id=qualifiers.pop('id', '<unknown id>'),
                              strand=int(strands[row['strand']]),
                              qualifiers=qualifiers))
        return feature_list
