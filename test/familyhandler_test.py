#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments
from __future__ import print_function, absolute_import, division
from argparse import Namespace

from connor.familyhandler import _CigarStatHandler, build_family_handlers
from connor.familyhandler import _CigarMinorityStatHandler
from connor.familyhandler import _FamilySizeStatHandler
from connor.familyhandler import _MatchStatHandler
from connor.familyhandler import _WriteAnnotatedAlignsHandler
from connor.familyhandler import _WriteConsensusHandler
import connor.samtools as samtools
from connor.samtools import ConnorAlign
from test.connor_test import _mock_tag_family
from test.samtools_test import mock_align
import test.samtools_test as samtools_test
from test.utils_test import BaseConnorTestCase
from test.utils_test import MicroMock

class FamilyHandlerTest(BaseConnorTestCase):
    def test_build_family_handlers(self):
        args = Namespace(umt_distance_threshold=1,
                         min_family_size_threshold=3)

        handlers = build_family_handlers(args,
                                         samtools_test.MockAlignWriter(),
                                         samtools.AlignWriter.NULL,
                                         self.mock_logger)
        actual_handler_names = [x.__class__.__name__ for x in handlers]

        expected_handler_names = ['_FamilySizeStatHandler',
                                  '_MatchStatHandler',
                                  '_CigarMinorityStatHandler',
                                  '_CigarStatHandler',
                                  '_WriteConsensusHandler']
        self.assertEqual(expected_handler_names, actual_handler_names)

    def test_build_family_handlers_withAnnotatedAlign(self):
        args = Namespace(umt_distance_threshold=1,
                         min_family_size_threshold=3)
        handlers = build_family_handlers(args,
                                         samtools_test.MockAlignWriter(),
                                         samtools_test.MockAlignWriter(),
                                         self.mock_logger)
        actual_handler_names = [x.__class__.__name__ for x in handlers]

        expected_handler_names = ['_FamilySizeStatHandler',
                                  '_MatchStatHandler',
                                  '_CigarMinorityStatHandler',
                                  '_CigarStatHandler',
                                  '_WriteConsensusHandler',
                                  '_WriteAnnotatedAlignsHandler']
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

    def test_end_median(self):
        posAfam1 = MicroMock(align_pairs=[1,1])
        posAfam2 = MicroMock(align_pairs=[1,1,1])
        posBfam1 = MicroMock(align_pairs=[1,1,1,1,1])
        families = [posAfam1, posAfam2, posBfam1]
        stat_handler = _FamilySizeStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(3, stat_handler.median)

    def test_end_quantiles(self):
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

    def test_summary(self):
        stat_handler = _FamilySizeStatHandler(self.mock_logger)
        stat_handler.min = 1
        stat_handler.quartile_1 = 2
        stat_handler.median = 3
        stat_handler.quartile_3 = 4
        stat_handler.max = 5

        self.assertEqual((1,2,3,4,5), stat_handler.summary)

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

def _mock_align_pair(query_name, filter_value=None):
    left = ConnorAlign(mock_align(query_name=query_name), filter_value)
    right = ConnorAlign(mock_align(query_name=query_name), filter_value)
    return MicroMock(left=left, right=right)

def _mock_align_pairs(num_pairs, query_prefix):
    pairs = []
    query_name_fmt = "{}_{}"
    for i in range(0, num_pairs):
        pairs.append(_mock_align_pair(query_name_fmt.format(query_prefix, i)))
    return pairs

class WriteConsensusHandlerTest(BaseConnorTestCase):
    def test_handle_writesConsensus(self):
        family_1 = _mock_tag_family(consensus=_mock_align_pair("readA"),
                                    filter_value=None)
        family_2 = _mock_tag_family(consensus=_mock_align_pair("readB"),
                                    filter_value=None)
        families = [family_1, family_2]
        writer = samtools_test.MockAlignWriter()

        handler = _WriteConsensusHandler(writer)
        for family in families:
            handler.handle(family)

        name_pairs = [(fam,
                       align.query_name) for fam, align in writer._write_calls]
        self.assertEqual([(family_1, 'readA'), (family_1, 'readA'),
                          (family_2, 'readB'), (family_2, 'readB')], name_pairs)

    def test_handle_excludeFilteredFamilies(self):
        family_1 = _mock_tag_family(consensus=_mock_align_pair("readA"),
                                    filter_value=None)
        family_2 = _mock_tag_family(consensus=_mock_align_pair("readB"),
                                    filter_value='foo')
        families = [family_1, family_2]
        writer = samtools_test.MockAlignWriter()

        handler = _WriteConsensusHandler(writer)
        for family in families:
            handler.handle(family)

        name_pairs = [(fam,
                       align.query_name) for fam, align in writer._write_calls]
        self.assertEqual([(family_1, 'readA'), (family_1, 'readA')], name_pairs)

class WriteAnnotatedAlignsHandlerTest(BaseConnorTestCase):
    def test_handle_writesAllAlignments(self):
        pairA1 = _mock_align_pair("readA1")
        pairA2 = _mock_align_pair("readA2", filter_value="foo")
        pairB1 = _mock_align_pair("readB1")
        pairB2 = _mock_align_pair("readB2", filter_value="bar")

        family_A = _mock_tag_family(align_pairs=[pairA1, pairA2])
        family_B = _mock_tag_family(align_pairs=[pairB1, pairB2])
        families = [family_A, family_B]
        writer = samtools_test.MockAlignWriter()

        handler = _WriteAnnotatedAlignsHandler(writer)
        for family in families:
            handler.handle(family)

        name_pairs = [(fam,
                       align.query_name) for fam, align in writer._write_calls]
        self.assertEqual([(family_A, 'readA1'), (family_A, 'readA1'),
                          (family_A, 'readA2'), (family_A, 'readA2'),
                          (family_B, 'readB1'), (family_B, 'readB1'),
                          (family_B, 'readB2'), (family_B, 'readB2')], name_pairs)

class CigarsStatHandlerTest(BaseConnorTestCase):
    def test_percent_deduplication(self):
        posAfam1 = _mock_tag_family(included_pair_count=1,
                                   align_pairs = [1] * 1)
        posAfam2 = _mock_tag_family(included_pair_count=2,
                                   align_pairs = [1] * 2)
        posAfam3 = _mock_tag_family(included_pair_count=2,
                                   align_pairs = [1] * 2)
        families = [posAfam1, posAfam2, posAfam3]
        stat_handler = _CigarStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        expected = 1 - (3/5)
        self.assertEqual(expected, stat_handler.percent_deduplication)

    def test_total_input_alignment_allCounted(self):
        posAfam1 = _mock_tag_family(included_pair_count=1,
                                   align_pairs = [1] * 1)
        posAfam2 = _mock_tag_family(included_pair_count=2,
                                   align_pairs = [1] * 2)
        families = [posAfam1, posAfam2]
        stat_handler = _CigarStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(3, stat_handler.total_input_alignment_count)
        self.assertEqual(3, stat_handler.total_alignment_count)
        self.assertEqual(0, stat_handler.total_excluded_alignments)
        self.assertEqual(0, stat_handler.percent_excluded_alignments)

    def test_total_input_alignment_someCounted(self):
        posAfam1 = _mock_tag_family(included_pair_count=3,
                                   align_pairs = [1] * 3)
        posAfam2 = _mock_tag_family(included_pair_count=1,
                                   align_pairs = [1] * 2)
        families = [posAfam1, posAfam2]
        stat_handler = _CigarStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(5, stat_handler.total_input_alignment_count)
        self.assertEqual(4, stat_handler.total_alignment_count)
        self.assertEqual(1, stat_handler.total_excluded_alignments)
        self.assertEqual(0.2, stat_handler.percent_excluded_alignments)

    def test_end_min(self):
        posAfam1 = _mock_tag_family(distinct_cigar_count=1)
        posAfam2 = _mock_tag_family(distinct_cigar_count=2)
        posBfam1 = _mock_tag_family(distinct_cigar_count=4)
        families = [posAfam1, posAfam2, posBfam1]
        stat_handler = _CigarStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(1, stat_handler.min)

    def test_end_max(self):
        posAfam1 = _mock_tag_family(distinct_cigar_count=1)
        posAfam2 = _mock_tag_family(distinct_cigar_count=2)
        posBfam1 = _mock_tag_family(distinct_cigar_count=4)
        families = [posAfam1, posAfam2, posBfam1]
        stat_handler = _CigarStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(4, stat_handler.max)

    def test_end_median(self):
        posAfam1 = _mock_tag_family(distinct_cigar_count=1)
        posAfam2 = _mock_tag_family(distinct_cigar_count=2)
        posBfam1 = _mock_tag_family(distinct_cigar_count=4)
        families = [posAfam1, posAfam2, posBfam1]
        stat_handler = _CigarStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(2, stat_handler.median)

    def test_end_quantiles(self):
        posAfam1 = _mock_tag_family(distinct_cigar_count=2)
        posAfam2 = _mock_tag_family(distinct_cigar_count=3)
        posBfam1 = _mock_tag_family(distinct_cigar_count=9)
        posBfam2 = _mock_tag_family(distinct_cigar_count=12)
        posBfam3 = _mock_tag_family(distinct_cigar_count=7)
        families = [posAfam1, posAfam2, posBfam1, posBfam2, posBfam3]
        stat_handler = _CigarStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(3, stat_handler.quartile_1)
        self.assertEqual(9, stat_handler.quartile_3)

    def test_summary(self):
        stat_handler = _CigarStatHandler(self.mock_logger)
        stat_handler.min = 1
        stat_handler.quartile_1 = 2
        stat_handler.median = 3
        stat_handler.quartile_3 = 4
        stat_handler.max = 5

        self.assertEqual((1,2,3,4,5), stat_handler.summary)


class CigarMinorityStatHandlerTest(BaseConnorTestCase):
    def test_minority_cigar_percentage_summary(self):
        posAfam0 = _mock_tag_family(minority_cigar_percentage=0)
        posAfam1 = _mock_tag_family(minority_cigar_percentage=0.01)
        posAfam2 = _mock_tag_family(minority_cigar_percentage=0.1)
        posAfam3 = _mock_tag_family(minority_cigar_percentage=0.3)
        posAfam4 = _mock_tag_family(minority_cigar_percentage=0.7)
        posAfam5 = _mock_tag_family(minority_cigar_percentage=0.9)
        families = [posAfam0,
                    posAfam1,
                    posAfam2,
                    posAfam3,
                    posAfam4,
                    posAfam5]

        stat_handler = _CigarMinorityStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual((0.01, 0.1, 0.3, 0.7, 0.9), stat_handler.summary)
