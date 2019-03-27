
#!/bin/env python
# -*- coding: utf-8 -*-
##
# test_serialization.py: Checks correctness of JSON serialization.
##
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License.
##

## IMPORTS ##

import unittest
import json
from qsharp.serialization import map_tuples, unmap_tuples

class TestSerialization(unittest.TestCase):
    def test_map_shallow_tuple(self):
        self.assertEqual(
            map_tuples((42, 'foo')),
            {'@type': 'tuple', 'item1': 42, 'item2': 'foo'}
        )

    def test_map_deep_tuple(self):
        actual = {
            'foo': [1, 3.14, (42, 'baz')],
            'bar': {'a': ('a', 'a'), 'b': ()}
        }
        expected = {
            'foo': [1, 3.14, {'@type': 'tuple', 'item1': 42, 'item2': 'baz'}],
            'bar': {
                'a': {'@type': 'tuple', 'item1': 'a', 'item2': 'a'},
                'b': {'@type': 'tuple'}
            }
        }
        self.assertEqual(
            map_tuples(actual), expected
        )

    def test_roundtrip_shallow_tuple(self):
        actual = ('a', 3.14, False)
        self.assertEqual(
            unmap_tuples(map_tuples(actual)), actual
        )

    def test_roundtrip_dict(self):
        actual = {'a': 'b', 'c': ('d', 'e')}
        self.assertEqual(
            unmap_tuples(map_tuples(actual)), actual
        )

    def test_roundtrip_deep_tuple(self):
        actual = ('a', ('b', 'c'))
        self.assertEqual(
            unmap_tuples(map_tuples(actual)), actual
        )

    def test_roundtrip_very_deep_tuple(self):
        actual = {
            'a': {
                'b': (
                    {
                        'c': ('d', ['e', ('f', 'g', 12, False)])
                    },
                    ['h', {'g': ('i', 'j')}]
                )
            }
        }
        self.assertEqual(
            unmap_tuples(map_tuples(actual)), actual
        )

if __name__ == "__main__":
    unittest.main()
