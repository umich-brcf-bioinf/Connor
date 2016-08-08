#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments
from __future__ import print_function, absolute_import, division
from argparse import Namespace
import os

from testfixtures.tempdirectory import TempDirectory

from test.utils_test import BaseConnorTestCase
from test.utils_test import MicroMock
from test.samtools_test import mock_align
from connor.familyhandler import _CigarStatHandler, build_family_handlers
from connor.familyhandler import _CigarMinorityStatHandler
from connor.familyhandler import _FamilySizeStatHandler
from connor.familyhandler import _MatchStatHandler
from connor.familyhandler import _WriteAnnotatedAlignsHandler
import connor.samtools as samtools
import test.samtools_test as samtools_test

def _mock_tag_family(input_alignment_count=5,
                    alignments=None,
                    excluded_alignments=None,
                    distinct_cigar_count=1,
                    inexact_match_count=0,
                    minority_cigar_percentage=0):
    if alignments is None:
        alignments = [1,2,3,4]
    if excluded_alignments is None:
        excluded_alignments = [5,6]
    return MicroMock(input_alignment_count=input_alignment_count,
                     alignments=alignments,
                     excluded_alignments=excluded_alignments,
                     distinct_cigar_count=distinct_cigar_count,
                     inexact_match_count=inexact_match_count,
                     minority_cigar_percentage=minority_cigar_percentage)


class FamilyHandlerTest(BaseConnorTestCase):
    def test_build_family_handlers(self):
        args = Namespace(umi_distance_threshold=1,
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
                                  '_WriteFamilyHandler']
        self.assertEqual(expected_handler_names, actual_handler_names)

    def test_build_family_handlers_withAnnotatedAlign(self):
        args = Namespace(umi_distance_threshold=1,
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
                                  '_WriteFamilyHandler',
                                  '_WriteAnnotatedAlignsHandler']
        self.assertEqual(expected_handler_names, actual_handler_names)


class FamilySizeStatHandlerTest(BaseConnorTestCase):
    def test_end_min(self):
        posAfam1 = MicroMock(alignments=[1,1])
        posAfam2 = MicroMock(alignments=[1,1,1])
        posBfam1 = MicroMock(alignments=[1,1,1,1,1])
        families = [posAfam1, posAfam2, posBfam1]
        stat_handler = _FamilySizeStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(2, stat_handler.min)

    def test_end_max(self):
        posAfam1 = MicroMock(alignments=[1,1])
        posAfam2 = MicroMock(alignments=[1,1,1])
        posBfam1 = MicroMock(alignments=[1,1,1,1,1])
        families = [posAfam1, posAfam2, posBfam1]
        stat_handler = _FamilySizeStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(5, stat_handler.max)

    def test_end_median(self):
        posAfam1 = MicroMock(alignments=[1,1])
        posAfam2 = MicroMock(alignments=[1,1,1])
        posBfam1 = MicroMock(alignments=[1,1,1,1,1])
        families = [posAfam1, posAfam2, posBfam1]
        stat_handler = _FamilySizeStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(3, stat_handler.median)

    def test_end_quantiles(self):
        posAfam1 = MicroMock(alignments=[1]*2)
        posAfam2 = MicroMock(alignments=[1]*3)
        posBfam1 = MicroMock(alignments=[1]*9)
        posBfam2 = MicroMock(alignments=[1]*12)
        posBfam3 = MicroMock(alignments=[1]*8)
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
        args = Namespace(umi_distance_threshold=1)
        stat_handler = _MatchStatHandler(args, self.mock_logger)
        posAfam1 = _mock_tag_family(alignments=[1] * 5,
                                   inexact_match_count=1)
        posAfam2 = _mock_tag_family(alignments=[1] * 15,
                                   inexact_match_count=4)
        families = [posAfam1, posAfam2]

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        self.assertEqual(5, stat_handler.total_inexact_match_count)
        self.assertEqual(20, stat_handler.total_pair_count)
        self.assertEqual(5/20, stat_handler.percent_inexact_match)

def align_pairs(num_pairs, query_prefix):
    pairs = []
    query_name_format = "{}_{}"
    for i in range(0, num_pairs):
        left = mock_align(query_name = query_name_format.format(query_prefix, i))
        right = mock_align(query_name = query_name_format.format(query_prefix, i))
        pair =  MicroMock(left_alignment=left, right_alignment=right)
        pairs.append(pair)
    return pairs

class WriteAnnotatedAlignsHandlerTest(BaseConnorTestCase):
    def test_handle(self):
        family_1 = _mock_tag_family(alignments=align_pairs(1, "A"),
                                    excluded_alignments=align_pairs(2, "a"))
        family_2 = _mock_tag_family(alignments=align_pairs(1, "B"),
                                    excluded_alignments=align_pairs(2, "b"))
        families = [family_1, family_2]
        writer = samtools_test.MockAlignWriter()

        handler = _WriteAnnotatedAlignsHandler(writer)
        for family in families:
            handler.handle(family)

        name_pairs = [(fam,
                       align.query_name) for fam, align in writer._write_calls]
        self.assertEqual([(family_1, 'A_0'), (family_1, 'A_0'),
                          (family_1, 'a_0'), (family_1, 'a_0'),
                          (family_1, 'a_1'), (family_1, 'a_1'),
                          (family_2, 'B_0'), (family_2, 'B_0'),
                          (family_2, 'b_0'), (family_2, 'b_0'),
                          (family_2, 'b_1'), (family_2, 'b_1')], name_pairs)

class CigarsStatHandlerTest(BaseConnorTestCase):
    def test_percent_deduplication(self):
        posAfam1 = _mock_tag_family(input_alignment_count=1,
                                   alignments = [1] * 1)
        posAfam2 = _mock_tag_family(input_alignment_count=2,
                                   alignments = [1] * 2)
        posAfam3 = _mock_tag_family(input_alignment_count=2,
                                   alignments = [1] * 2)
        families = [posAfam1, posAfam2, posAfam3]
        stat_handler = _CigarStatHandler(self.mock_logger)

        for family in families:
            stat_handler.handle(family)
        stat_handler.end()

        expected = 1 - (3/5)
        self.assertEqual(expected, stat_handler.percent_deduplication)

    def test_total_input_alignment_allCounted(self):
        posAfam1 = _mock_tag_family(input_alignment_count=1,
                                   alignments = [1] * 1)
        posAfam2 = _mock_tag_family(input_alignment_count=2,
                                   alignments = [1] * 2)
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
        posAfam1 = _mock_tag_family(input_alignment_count=3,
                                   alignments = [1] * 3)
        posAfam2 = _mock_tag_family(input_alignment_count=2,
                                   alignments = [1] * 1)
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
