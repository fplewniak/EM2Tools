import argparse


class get_list(argparse.Action):
    def __init__(self, option_strings, dest, nargs='+', **kwargs):
        if nargs is not '+':
            raise ValueError('get_list action requires nargs to be set to "+"')
        super(get_list, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, self.arg2list(values))


def arg2list(values):
    arglist = []
    if values is not None:
        for v in values:
            try:
                listfile = open(v)
                arglist += [line.strip('\n') for line in listfile]
                listfile.close()
            except IOError:
                arglist.append(v)
    return arglist
