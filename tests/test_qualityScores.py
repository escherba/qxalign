#!/usr/bin/env python3

import unittest
from qxalign import Qxalign


class TestQualityScores(unittest.TestCase):

    def test_qualityScores(self):
        q = Qxalign()
        q.prepare("AAAACGT", "TGCA", "!!!!!!!!!!!")
        self.assertEqual(60, q.align())

    def test_qualityScoresTrace(self):
        q = Qxalign()

        q.prepare("AAAACGT", "TGCA", b"!!!!")
        self.assertEqual(60, q.align())
        q.trace()
        self.assertEqual("3I 1=", q.show_trace())

        q.prepare_query(query_seq="CAAC")
        self.assertEqual(40, q.align(semi=True))
        q.trace()
        self.assertEqual("1X 3=", q.show_trace())

        q.prepare("", "", "")
        self.assertRaises(IndexError, q.align, [])


if __name__ == "__main__":
    unittest.run(verbose=True)
