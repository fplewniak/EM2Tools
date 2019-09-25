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
        elif code.upper() in 'BZX':
            return {'B': 'DN', 'Z': 'EQ', 'X': '.'}[code.upper()]
        else:
            raise ValueError('Invalid amino-acid %s' % code)
    elif code.upper() in IUPACData.ambiguous_dna_values:
        return ''.join(sorted(IUPACData.ambiguous_dna_values[code.upper()]))
    else:
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
    pattern = re.sub(r'{(\w+)}', r'[^\1]', pattern)           # none of the specified residues
    pattern = re.sub(r'^<', r'^', pattern)                    # start of sequence
    pattern = re.sub(r'>$', r'$', pattern)                    # end of sequence
    pattern = re.sub(r'\((\d+(,\d*)*)\)', r'{\1}', pattern)   # repeat residue or subsequence

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
