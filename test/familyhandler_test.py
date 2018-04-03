#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments
from __future__ import print_function, absolute_import, division
from argparse import Namespace

from test.connor_test import _mock_tag_family
from test.consam_test.writers_test import MockAlignWriter
from test.utils_test import BaseConnorTestCase
from test.utils_test import MicroMock

from connor.familyhandler import build_family_handlers
from connor.familyhandler import _FamilySizeStatHandler
from connor.familyhandler import _MatchStatHandler
from connor.familyhandler import _WriteAnnotatedAlignsHandler
from connor.familyhandler import _WriteConsensusHandler
from connor.consam.alignments import ConnorAlign


class FamilyHandlerTest(BaseConnorTestCase):
    def test_build_family_handlers(self):
        args = Namespace(umt_distance_threshold=1,
                         min_family_size_threshold=3)
        handlers = build_family_handlers(args,
                                         MockAlignWriter(),
                                         MockAlignWriter(),
                                         self.mock_logger)
        actual_handler_names = [x.__class__.__name__ for x in handlers]

        expected_handler_names = ['_FamilySizeStatHandler',
                                  '_MatchStatHandler',
                                  '_WriteAnnotatedAlignsHandler',
                                  '_WriteConsensusHandler']
        self.assertEqual(expected_handler_names, actual_handler_names)


class FamilySizeStatHandlerTest(BaseConnorTestCase):
    def test_end_min(self):
        posAfam1 = MicroMock(align_pairs=[1,1])
        posAfam2 = MicroMock(align_pairs=[1,1,1])
        posBfam1 = MicroMock(align_pairs=[1,1,1,1,1])
        families = [posAfam1, posAfam2, posBfam1]
        stat_handler = _FamilySizeStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(2, stat_handler.min)

    def test_end_max(self):
        posAfam1 = MicroMock(align_pairs=[1,1])
        posAfam2 = MicroMock(align_pairs=[1,1,1])
        posBfam1 = MicroMock(align_pairs=[1,1,1,1,1])
        families = [posAfam1, posAfam2, posBfam1]
        stat_handler = _FamilySizeStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(5, stat_handler.max)

    def test_end_median_odd(self):
        posAfam1 = MicroMock(align_pairs=[1,1])
        posAfam2 = MicroMock(align_pairs=[1,1,1])
        posBfam1 = MicroMock(align_pairs=[1,1,1,1,1])
        families = [posAfam1, posAfam2, posBfam1]
        stat_handler = _FamilySizeStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(3, stat_handler.median)

    def test_end_median_even(self):
        posAfam1 = MicroMock(align_pairs=[1,1])
        posAfam2 = MicroMock(align_pairs=[1,1,1])
        posBfam1 = MicroMock(align_pairs=[1,1,1,1])
        posBfam2 = MicroMock(align_pairs=[1,1,1,1,1])
        families = [posAfam1, posAfam2, posBfam1, posBfam2]
        stat_handler = _FamilySizeStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(3.5, stat_handler.median)

    def test_end_mean(self):
        posAfam1 = MicroMock(align_pairs=[1] * 1)
        posAfam2 = MicroMock(align_pairs=[1] * 2)
        posBfam1 = MicroMock(align_pairs=[1] * 4)
        posBfam2 = MicroMock(align_pairs=[1] * 8)
        posBfam3 = MicroMock(align_pairs=[1] * 16)
        families = [posAfam1, posAfam2, posBfam1, posBfam2, posBfam3]
        stat_handler = _FamilySizeStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(6.2, stat_handler.mean)

    def test_end_quantiles_odd(self):
        posAfam1 = MicroMock(align_pairs=[1]*2)
        posAfam2 = MicroMock(align_pairs=[1]*3)
        posBfam1 = MicroMock(align_pairs=[1]*9)
        posBfam2 = MicroMock(align_pairs=[1]*12)
        posBfam3 = MicroMock(align_pairs=[1]*8)
        families = [posAfam1, posAfam2, posBfam1, posBfam2, posBfam3]
        stat_handler = _FamilySizeStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(3, stat_handler.quartile_1)
        self.assertEqual(9, stat_handler.quartile_3)

    def test_end_quantiles_even(self):
        posAfam1 = MicroMock(align_pairs=[1]*2)
        posAfam2 = MicroMock(align_pairs=[1]*3)
        posBfam1 = MicroMock(align_pairs=[1]*9)
        posBfam2 = MicroMock(align_pairs=[1]*12)
        families = [posAfam1, posAfam2, posBfam1, posBfam2]
        stat_handler = _FamilySizeStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(2.75, stat_handler.quartile_1)
        self.assertEqual(9.75, stat_handler.quartile_3)

    def test_summary(self):
        stat_handler = _FamilySizeStatHandler(self.mock_logger)
        stat_handler.min = 1
        stat_handler.quartile_1 = 2
        stat_handler.median = 4
        stat_handler.mean = 6.2
        stat_handler.quartile_3 = 8
        stat_handler.max = 16

        self.assertEqual((1, 2, 4, 6.2, 8, 16), stat_handler.summary)

    def test_end_zeroFamilies(self):
        stat_handler = _FamilySizeStatHandler(self.mock_logger)
        stat_handler.end()

        self.assertEqual(None, stat_handler.min)
        self.assertEqual(None, stat_handler.quartile_1)
        self.assertEqual(None, stat_handler.median)
        self.assertEqual(None, stat_handler.mean)
        self.assertEqual(None, stat_handler.quartile_3)
        self.assertEqual(None, stat_handler.max)

class MatchStatHandlerTest(BaseConnorTestCase):
    def test_total_inexact_match_count(self):
        args = Namespace(umt_distance_threshold=1)
        stat_handler = _MatchStatHandler(args, self.mock_logger)
        posAfam1 = _mock_tag_family(align_pairs=[1] * 5,
                                   inexact_match_count=1)
        posAfam2 = _mock_tag_family(align_pairs=[1] * 15,
                                   inexact_match_count=4)
        families = [posAfam1, posAfam2]

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(5, stat_handler.total_inexact_match_count)
        self.assertEqual(20, stat_handler.total_pair_count)
        self.assertEqual(5/20, stat_handler.percent_inexact_match)

    def test_total_percent_inexact_match_whenZeroFamilies(self):
        args = Namespace(umt_distance_threshold=1)
        stat_handler = _MatchStatHandler(args, self.mock_logger)
        stat_handler.end()

        self.assertEqual(0, stat_handler.total_inexact_match_count)
        self.assertEqual(0, stat_handler.total_pair_count)
        self.assertEqual(0, stat_handler.percent_inexact_match)


class WriteConsensusHandlerTest(BaseConnorTestCase):
    def _mock_align_pair(self, query_name, filter_value=None):
        left = ConnorAlign(self.mock_align(query_name=query_name), filter_value)
        right = ConnorAlign(self.mock_align(query_name=query_name), filter_value)
        return MicroMock(query_name=query_name, left=left, right=right)

    def _mock_align_pairs(self, num_pairs, query_prefix):
        pairs = []
        query_name_fmt = "{}_{}"
        for i in range(0, num_pairs):
            pairs.append(self._mock_align_pair(query_name_fmt.format(query_prefix, i)))
        return pairs

    def test_handle_writesConsensus(self):
        family_1 = _mock_tag_family(consensus=self._mock_align_pair("readA"),
                                    filter_value=None)
        family_2 = _mock_tag_family(consensus=self._mock_align_pair("readB"),
                                    filter_value=None)
        families = [family_1, family_2]
        writer = MockAlignWriter()

        handler = _WriteConsensusHandler(writer)
        for family in families:
            handler.handle(family)

        name_pairs = [(fam,
                       align.query_name) for fam, align in writer._write_calls]
        self.assertEqual([(family_1, 'readA'), (family_1, 'readA'),
                          (family_2, 'readB'), (family_2, 'readB')], name_pairs)

    def test_handle_excludeFilteredFamilies(self):
        family_1 = _mock_tag_family(consensus=self._mock_align_pair("readA"),
                                    filter_value=None)
        family_2 = _mock_tag_family(consensus=self._mock_align_pair("readB"),
                                    filter_value='foo')
        families = [family_1, family_2]
        writer = MockAlignWriter()

        handler = _WriteConsensusHandler(writer)
        for family in families:
            handler.handle(family)

        name_pairs = [(fam,
                       align.query_name) for fam, align in writer._write_calls]
        self.assertEqual([(family_1, 'readA'), (family_1, 'readA')], name_pairs)

class WriteAnnotatedAlignsHandlerTest(BaseConnorTestCase):
    def _mock_align_pair(self, query_name, filter_value=None):
        left = ConnorAlign(self.mock_align(query_name=query_name), filter_value)
        right = ConnorAlign(self.mock_align(query_name=query_name), filter_value)
        return MicroMock(query_name=query_name, left=left, right=right)

    def test_handle_writesAllAlignments(self):
        pairA1 = self._mock_align_pair("readA1")
        pairA2 = self._mock_align_pair("readA2", filter_value="foo")
        pairB1 = self._mock_align_pair("readB1")
        pairB2 = self._mock_align_pair("readB2", filter_value="bar")

        family_A = _mock_tag_family(align_pairs=set([pairA1, pairA2]))
        family_B = _mock_tag_family(align_pairs=set([pairB1, pairB2]))
        families = [family_A, family_B]
        writer = MockAlignWriter()

        handler = _WriteAnnotatedAlignsHandler(writer)
        for family in families:
            handler.handle(family)

        name_pairs = [(fam,
                       align.query_name) for fam, align in writer._write_calls]
        self.assertEqual([(family_A, 'readA1'), (family_A, 'readA1'),
                          (family_A, 'readA2'), (family_A, 'readA2'),
                          (family_B, 'readB1'), (family_B, 'readB1'),
                          (family_B, 'readB2'), (family_B, 'readB2')],
                         name_pairs)

    def test_handle_ordersAlignmentsByReadname(self):
        pairA1 = self._mock_align_pair("readA1")
        pairA2 = self._mock_align_pair("readA2")
        pairA3 = self._mock_align_pair("readA3")
        pairA4 = self._mock_align_pair("readA4")

        family_A = _mock_tag_family(align_pairs=set([pairA2,
                                                     pairA1,
                                                     pairA4,
                                                     pairA3]))
        families = [family_A]
        writer = MockAlignWriter()

        handler = _WriteAnnotatedAlignsHandler(writer)
        for family in families:
            handler.handle(family)

        names = [align.query_name for _, align in writer._write_calls]
        self.assertEqual(['readA1', 'readA1',
                          'readA2', 'readA2',
                          'readA3', 'readA3',
                          'readA4', 'readA4'], names)
