"""
Tests for the extension module to the standard argparse module.
"""
import unittest
from em2lib import argparse_em2


class GetListAction(unittest.TestCase):
    # Test the GetList argparse custom GetList that returns a list from either a list of arguments and/or
    # file contents
    @staticmethod
    def test_fromfile():
        assert argparse_em2.GetList.arg2list(['data/list_1columnA.txt']) == ['Item1', 'Item2', 'Item3']

    @staticmethod
    def test_from2files():
        assert argparse_em2.GetList.arg2list(
            ['data/list_1columnA.txt', 'data/list_1columnB.txt']) == ['Item1', 'Item2', 'Item3', 'Item4', 'Item5']

    @staticmethod
    def test_from_one():
        assert argparse_em2.GetList.arg2list(['a']) == ['a']

    @staticmethod
    def test_fromlist():
        assert argparse_em2.GetList.arg2list(['a', 'b', 'c']) == ['a', 'b', 'c']

    @staticmethod
    def test_fromfile_and_list():
        assert argparse_em2.GetList.arg2list(['data/list_1columnA.txt', 'a', 'b', 'c']) == ['Item1', 'Item2', 'Item3',
                                                                                           'a', 'b', 'c']

    @staticmethod
    def test_none():
        assert argparse_em2.GetList.arg2list(None) == []


class GetSetAction(unittest.TestCase):
    # Test the GetSet argparse custom GetSet action that returns a set from either a list of arguments and/or
    # file contents
    @staticmethod
    def test_fromfile():
        assert argparse_em2.GetSet.arg2set(['data/list_1columnA.txt']) == {'Item1', 'Item2', 'Item3'}

    @staticmethod
    def test_from2files():
        assert argparse_em2.GetSet.arg2set(
            ['data/list_1columnA.txt', 'data/list_1columnB.txt']) == {'Item1', 'Item2', 'Item3', 'Item4', 'Item5'}

    @staticmethod
    def test_from2files_redundant():
        assert argparse_em2.GetSet.arg2set(
            ['data/list_1columnA.txt', 'data/list_1columnC.txt']) == {'Item1', 'Item2', 'Item3', 'Item4', 'Item5'}

    @staticmethod
    def test_from_one():
        assert argparse_em2.GetSet.arg2set(['a']) == {'a'}

    @staticmethod
    def test_fromlist():
        assert argparse_em2.GetSet.arg2set(['a', 'b', 'c']) == {'a', 'b', 'c'}

    @staticmethod
    def test_fromlist_redundant():
        assert argparse_em2.GetSet.arg2set(['a', 'b', 'c', 'a']) == {'a', 'b', 'c'}

    @staticmethod
    def test_fromfile_and_list():
        assert argparse_em2.GetSet.arg2set(['data/list_1columnA.txt', 'a', 'b', 'c']) == {'Item1', 'Item2', 'Item3',
                                                                                         'a', 'b', 'c'}

    @staticmethod
    def test_fromfile_and_list_redundant():
        assert argparse_em2.GetSet.arg2set(['data/list_1columnA.txt', 'a', 'Item2', 'c']) == {'Item1', 'Item2', 'Item3',
                                                                                             'a', 'c'}

    @staticmethod
    def test_none():
        assert argparse_em2.GetSet.arg2set(None) == set()
