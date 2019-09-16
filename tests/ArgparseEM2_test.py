'''
Tests for the extension module to the standard argparse module.
'''
import unittest
from EM2Libs import ArgparseEM2


class GetListAction(unittest.TestCase):
    # Test the GetList argparse custom get_list_action that returns a list from either a list of arguments and/or
    # file contents
    @staticmethod
    def test_fromfile():
        assert ArgparseEM2.GetList.arg2list(['data/list_1columnA.txt']) == ['Item1', 'Item2', 'Item3']

    @staticmethod
    def test_from2files():
        assert ArgparseEM2.GetList.arg2list(
            ['data/list_1columnA.txt', 'data/list_1columnB.txt']) == ['Item1', 'Item2', 'Item3', 'Item4', 'Item5']

    @staticmethod
    def test_from_one():
        assert ArgparseEM2.GetList.arg2list(['a']) == ['a']

    @staticmethod
    def test_fromlist():
        assert ArgparseEM2.GetList.arg2list(['a', 'b', 'c']) == ['a', 'b', 'c']

    @staticmethod
    def test_fromfile_and_list():
        assert ArgparseEM2.GetList.arg2list(['data/list_1columnA.txt', 'a', 'b', 'c']) == ['Item1', 'Item2', 'Item3',
                                                                                           'a', 'b', 'c']

    @staticmethod
    def test_none():
        assert ArgparseEM2.GetList.arg2list(None) == []
