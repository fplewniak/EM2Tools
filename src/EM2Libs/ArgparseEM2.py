"""
Extension module to the standard argparse module. Adding custom actions and argument verification methods.
"""
import argparse


class GetList(argparse.Action):
    """
An argparse custom action to return a list from an argument containing a list of elements and/or file names. Files are
supposed to contain one element of the list per line. There can be more than one file and the argument may take a
combination of elements and files. In all cases, the returned list will contain all the specified elements without any
checking for redundancy. If you need a non redundant set instead of a list, then use GetSet action instead.
    """

    def __init__(self, option_strings, dest, nargs='+', **kwargs):
        if nargs != '+':
            raise ValueError('GetList get_list_action requires nargs to be set to "+"')
        super(GetList, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, self.arg2list(values))

    @staticmethod
    def arg2list(values):
        """
This method converts the argument values containing elements and/or files containing elements into a list of elements.
        :param values: argument values, this is supposed to be a list of arguments or None (returns an empty list)
        :return: the list of elements or an empty list if the argument was None
        """
        arglist = []
        if values is not None:
            for v in values:
                try:
                    arglist += [line.strip('\n') for line in open(v)]
                except IOError:
                    arglist.append(v)
        return arglist
