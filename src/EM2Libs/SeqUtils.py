"""
Extension module to the Biopython Bio.SeqUtils module
"""
import re

from Bio.Data import IUPACData


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


def seqfilter(records, minlength=None, maxlength=None, pattern=None, name=None, keep=True):
    """
    Filters a list of SeqRecords instances, keeping only records satisfying the specified criteria of length, match
    of a pattern, name specification. It is possible to invert the filtering process by setting the keep boolean to
    False and thus only keep records which do not satisfy the criteria.
    :param records: list of SeqRecord instances to seqfilter
    :param minlength: minimum length of sequence
    :param maxlength: maximum length of sequence
    :param pattern: sequence pattern
    :param name: sequence name
    :param keep: boolean, if True, keep the records respecting the criteria, otherwise, discard them and keep the others
    """
    filtered = []
    for r in records:
        if keep is True:
            if r.seq.length_in_range(minlength, maxlength):
                if (pattern is None) or (r.seq.search(pattern) != []):
                    if (name is None) or (re.match(name, r.name) is not None):
                        filtered.append(r)
        else:
            if (minlength is None) and (maxlength is None) or r.seq.length_in_range(minlength, maxlength) is False:
                if (pattern is None) or (r.seq.search(pattern) == []):
                    if (name is None) or (re.match(name, r.name) is None):
                        filtered.append(r)

    return filtered
