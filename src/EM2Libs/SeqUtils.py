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
        else:
            return {'B': 'DN', 'Z': 'EQ', 'X': '.'}[code.upper()]
    else:
        return ''.join(sorted(IUPACData.ambiguous_dna_values[code.upper()]))


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
    pattern = re.sub(r'{(\w+)}',r'[^\1]', pattern)
    pattern = re.sub(r'^<',r'^', pattern)
    pattern = re.sub(r'>$',r'$', pattern)
    pattern = re.sub(r'\((\d+(,\d+)*)\)',r'{\1}', pattern)

    return pattern
