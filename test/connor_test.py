#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments
#pylint: disable=deprecated-method

from __future__ import print_function, absolute_import, division
from argparse import Namespace
from collections import namedtuple
import os
import sys
from testfixtures.tempdirectory import TempDirectory
import connor.connor as connor
from connor import samtools
import connor.utils as utils
import test.samtools_test as samtools_test
from test.utils_test import BaseConnorTestCase
from test.utils_test import MicroMock

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

class MockAlignSegment(object):
    #pylint: disable=too-many-instance-attributes
    def __init__(self,
                 query_name,
                 reference_name,
                 reference_start,
                 next_reference_start,
                 query_sequence='AAACCC',
                 query_qualities=None,
                 cigarstring='6M',
                 reference_end=None):
        if not reference_end:
            reference_end = reference_start + len(query_sequence)
        self.query_name = query_name
        self.reference_name = reference_name
        self.reference_start = reference_start
        self.next_reference_start = next_reference_start
        self.query_sequence = query_sequence
        if query_qualities is None:
            self.query_qualities = [25,25,25,25,25,25]
        else:
            self.query_qualities = query_qualities
        self.cigarstring = cigarstring
        self.reference_end = reference_end
        self.filter = None

    def __hash__(self):
        return hash(self.query_name)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def set_tag(self, name, value, tag_type):
        pass


def align_seg(query_name, #pylint: disable=dangerous-default-value
              reference_name,
              reference_start,
              next_reference_start,
              query_sequence='AAACCC',
              query_qualities=[25,25,25,25,25,25],
              cigarstring='6M',
              reference_end=None):
    return MockAlignSegment(query_name,
                            reference_name,
                            reference_start,
                            next_reference_start,
                            query_sequence,
                            query_qualities,
                            cigarstring,
                            reference_end)

def align_pair(q, rn, rs, nrs, s1, s2, tag_length=3):
    alignL = align_seg(q, rn, rs, nrs, s1)
    alignR = align_seg(q, rn, rs, nrs, s2)
    return connor._PairedAlignment(alignL, alignR, tag_length)


class PairedAlignmentTest(BaseConnorTestCase):
    def test_init(self):
        left_align = MockAlignSegment("alignA", 'chr1', 10, 20,
                                      "AAATTT" "GGGG", reference_end=14)
        right_align = MockAlignSegment("alignA", 'chr1', 20, 10,
                                       "TTTT" "CCCGGG" ,reference_end=24)
        tag_length = 6
        actual_paired_alignment = connor._PairedAlignment(left_align,
                                                         right_align,
                                                         tag_length)

        self.assertIs(left_align, actual_paired_alignment.left_alignment)
        self.assertIs(right_align, actual_paired_alignment.right_alignment)
        self.assertEquals(("AAATTT","CCCGGG"), actual_paired_alignment.umi)

    def test_eq(self):
        left = align_seg("alignA", 'chr1', 100, 200, "AAANNNNNNN")
        right = align_seg("alignA", 'chr1', 200, 100, "CCCNNNNNNN")
        other = align_seg("alignB", 'chr1', 100, 200, "AAANNNNNNN")

        base = connor._PairedAlignment(left, right)
        self.assertEquals(base, connor._PairedAlignment(left, right))
        self.assertNotEquals(base, connor._PairedAlignment(other, right))
        self.assertNotEquals(base, connor._PairedAlignment(left, other))

    def test_hash(self):
        left_A = MockAlignSegment("alignA", 'chr1', 100, 200, "AAATTT" "NNNN")
        right_A = MockAlignSegment("alignA", 'chr1', 200, 100, "NNNN" "CCCGGG")
        left_B = MockAlignSegment("alignA", 'chr1', 100, 200, "AAATTT" "NNNN")
        right_B = MockAlignSegment("alignA", 'chr1', 200, 100, "NNNN" "CCCGGG")

        actual_set = set()
        base = connor._PairedAlignment(left_A, right_A)
        actual_set.add(base)
        self.assertEquals(1, len(actual_set))

        actual_set.add(base)
        self.assertEquals(1, len(actual_set))

        actual_set.add(connor._PairedAlignment(left_A, right_A))
        self.assertEquals(1, len(actual_set))

        equivalent_pair = connor._PairedAlignment(left_B, right_B)
        actual_set.add(equivalent_pair)
        self.assertEquals(1, len(actual_set))

    def test_replace_umi(self):
        left_A = align_seg('alignA', 'chr1', 100, 200, 'AAAA' 'NNNN')
        right_A = align_seg('alignA', 'chr1', 200, 100, 'NNNN' 'CCCC')
        paired_align = connor._PairedAlignment(left_A, right_A, tag_length=4)

        paired_align.replace_umi(('GGGG','TTTT'))

        self.assertEquals('GGGG' 'NNNN',
                          paired_align.left_alignment.query_sequence)
        self.assertEquals('NNNN' 'TTTT',
                          paired_align.right_alignment.query_sequence)


class TagFamiliyTest(BaseConnorTestCase):
    def test_init(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'GGGNNN', 'TTTNNN')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        alignments = [pair1, pair2, pair3]
        inexact_match_count = 42

        input_umi = "AAANNNCCCNNN"
        actual_tag_family = connor._TagFamily(input_umi,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(input_umi, actual_tag_family.umi)
        self.assertEquals(alignments, actual_tag_family.alignments)
        self.assertEquals(42, actual_tag_family.inexact_match_count)


    def test_consensus_rewrites_umi(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'GGGNNN', 'NNNTTT')
        input_umis = ("AAA", "CCC")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umis,
                                             [pair1],
                                             inexact_match_count,
                                             consensus_threshold=0.6)
        actual_consensus = actual_tag_family.consensus

        self.assertEquals(input_umis, actual_consensus.umi)
        self.assertEquals('AAANNN',
                          actual_consensus.left_alignment.query_sequence)
        self.assertEquals('NNNCCC',
                          actual_consensus.right_alignment.query_sequence)

    def test_consensus_sequence_trivial_noop(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umis,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)
        actual_consensus_seq = actual_tag_family.consensus

        self.assertEquals("nnnGGG",
                          actual_consensus_seq.left_alignment.query_sequence)
        self.assertEquals("TTTnnn",
                          actual_consensus_seq.right_alignment.query_sequence)

    def test_consensus_sequence_majority_wins(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'TCTnnn')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'TCTnnn')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")

        actual_tag_family = connor._TagFamily(input_umis,
                                              alignments,
                                              inexact_match_count=1,
                                              consensus_threshold=0.5)
        consensus_pair = actual_tag_family.consensus

        self.assertEquals("nnnGTG",
                          consensus_pair.left_alignment.query_sequence)
        self.assertEquals("TCTnnn",
                          consensus_pair.right_alignment.query_sequence)

    def test_consensus_sequence_below_threshold_Ns(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'TCTnnn')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'TCTnnn')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        threshold = 0.8
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umis,
                                             alignments,
                                             inexact_match_count,
                                             threshold)
        consensus_pair = actual_tag_family.consensus

        self.assertEquals("nnnGNG",
                          consensus_pair.left_alignment.query_sequence)
        self.assertEquals("TNTnnn",
                          consensus_pair.right_alignment.query_sequence)

    def test_consensus_qualities_majority_vote(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left_alignment.query_qualities = [30, 30, 30, 30, 30, 30]
        pair1.right_alignment.query_qualities = [25, 25, 25, 25, 25, 25]
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left_alignment.query_qualities = [30, 20, 30, 30, 30, 30]
        pair2.right_alignment.query_qualities = [25, 15, 25, 15, 15, 15]
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left_alignment.query_qualities = [10, 20, 10, 20, 20 ,20]
        pair3.right_alignment.query_qualities = [5, 15, 5 , 15, 15, 15]
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        inexact_match_count = 0


        actual_tag_family = connor._TagFamily(input_umis,
                                               alignments,
                                               inexact_match_count,
                                               consensus_threshold=0.6)
        paired_consensus = actual_tag_family.consensus

        self.assertEquals([30, 20, 30, 30, 30, 30],
                          paired_consensus.left_alignment.query_qualities)
        self.assertEquals([25, 15, 25, 15, 15, 15],
                          paired_consensus.right_alignment.query_qualities)

    def test_consensus_uniform_cigars_admitted(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left_alignment.cigarstring = "3S3M"
        pair1.right_alignment.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left_alignment.cigarstring = "3S3M"
        pair2.right_alignment.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left_alignment.cigarstring = "3S3M"
        pair3.right_alignment.cigarstring = "3S3M"
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umis,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        included_names = [a.left_alignment.query_name for a in actual_tag_family.alignments]
        self.assertEquals(['alignA','alignB','alignC'], included_names)

        excluded_names = [a.left_alignment.query_name for a in actual_tag_family.excluded_alignments]
        self.assertEquals([], excluded_names)


    def test_consensus_minority_cigars_excluded(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left_alignment.cigarstring = "3S3M"
        pair1.right_alignment.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left_alignment.cigarstring = "3S3M"
        pair2.right_alignment.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left_alignment.cigarstring = "3S3I"
        pair3.right_alignment.cigarstring = "3S3M"
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umis,
                                              alignments,
                                              inexact_match_count,
                                              consensus_threshold=0.6)

        included_names = [a.left_alignment.query_name for a in actual_tag_family.alignments]
        self.assertEquals(['alignA','alignB'], included_names)

        excluded_names = [a.left_alignment.query_name for a in actual_tag_family.excluded_alignments]
        self.assertEquals(['alignC'], excluded_names)


    def test_init_distinctCigarCountTwo(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left_alignment.cigarstring = "3S3M"
        pair1.right_alignment.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left_alignment.cigarstring = "3S3M"
        pair2.right_alignment.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left_alignment.cigarstring = "3S3I"
        pair3.right_alignment.cigarstring = "3S3I"
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umis,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(2, actual_tag_family.distinct_cigar_count)

    def test_init_distinctCigarCountOneForConsisentCigar(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left_alignment.cigarstring = "3S3M"
        pair1.right_alignment.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left_alignment.cigarstring = "3S3M"
        pair2.right_alignment.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left_alignment.cigarstring = "3S3M"
        pair3.right_alignment.cigarstring = "3S3M"
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umis,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(1, actual_tag_family.distinct_cigar_count)

    def test_distinct_cigar_count_single_ended(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left_alignment.cigarstring = "3S3M"
        pair1.right_alignment.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left_alignment.cigarstring = "3S3M"
        pair2.right_alignment.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left_alignment.cigarstring = "3S3I"
        pair3.right_alignment.cigarstring = "3S3M"
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umis,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(2, actual_tag_family.distinct_cigar_count)

    def test_init_noMinorityIfConsistentCigar(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left_alignment.cigarstring = "3S3M"
        pair1.right_alignment.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left_alignment.cigarstring = "3S3M"
        pair2.right_alignment.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left_alignment.cigarstring = "3S3M"
        pair3.right_alignment.cigarstring = "3S3M"
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umis,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(0, actual_tag_family.minority_cigar_percentage)

    def test_init_setsMinorityPercentage(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left_alignment.cigarstring = "3S3M"
        pair1.right_alignment.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left_alignment.cigarstring = "3S3M"
        pair2.right_alignment.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left_alignment.cigarstring = "3S3I"
        pair3.right_alignment.cigarstring = "3S3I"
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umis,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(1/3, actual_tag_family.minority_cigar_percentage)

    def test_init_addsFilterForMinority(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left_alignment.cigarstring = "3S3M"
        pair1.right_alignment.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left_alignment.cigarstring = "3S3M"
        pair2.right_alignment.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left_alignment.cigarstring = "3S3I"
        pair3.right_alignment.cigarstring = "3S3I"
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umis,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6,
                                             min_family_size=1)

        for a_pair in actual_tag_family.alignments:
            self.assertEqual(None, a_pair.left_alignment.filter)
            self.assertEqual(None, a_pair.right_alignment.filter)

        self.assertEquals(1, len(actual_tag_family.excluded_alignments))
        excluded_align_pair = actual_tag_family.excluded_alignments[0]
        self.assertEquals("minority CIGAR",
                          excluded_align_pair.left_alignment.filter)
        self.assertEquals("minority CIGAR",
                          excluded_align_pair.right_alignment.filter)

    def test_init_addsFilterForSmallFamily(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left_alignment.cigarstring = "3S3M"
        pair1.right_alignment.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left_alignment.cigarstring = "3S3M"
        pair2.right_alignment.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left_alignment.cigarstring = "3S3M"
        pair3.right_alignment.cigarstring = "3S3M"
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umis,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6,
                                             min_family_size=5)

        for a_pair in actual_tag_family.alignments:
            self.assertEqual('small family (<5)', a_pair.left_alignment.filter)
            self.assertEqual('small family (<5)', a_pair.right_alignment.filter)


    def test_input_alignment_count(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left_alignment.cigarstring = "3S3M"
        pair1.right_alignment.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left_alignment.cigarstring = "3S3M"
        pair2.right_alignment.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left_alignment.cigarstring = "3S3I"
        pair3.right_alignment.cigarstring = "3S3M"
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umis,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(3, actual_tag_family.input_alignment_count)


class ConnorTest(BaseConnorTestCase):
    @staticmethod
    def get_tag(tags, name):
        for tag in tags:
            if tag._tag_name == name:
                return tag
        return None

    def test_log_environment(self):
        args = Namespace(original_command_line=['foo', 'bar'],
                         other_stuff='baz')
        connor._log_environment_info(self.mock_logger, args)
        log_text = '\n'.join(self.mock_logger._log_calls['DEBUG'])
        self.assertRegexpMatches(log_text, 'command_line.*foo bar')
        self.assertRegexpMatches(log_text, 'command_options.*foo.*bar')
        self.assertRegexpMatches(log_text, 'command_options.*baz')
        self.assertRegexpMatches(log_text, 'command_cwd')
        self.assertRegexpMatches(log_text, 'platform_uname')
        self.assertRegexpMatches(log_text, 'python_version')
        self.assertRegexpMatches(log_text, 'pysam_version')

    def test_build_bam_tags(self):
        actual_tags = connor._build_bam_tags()
        self.assertEqual(5, len(actual_tags))

    def test_build_bam_tags_x0_filter(self):
        tag = ConnorTest.get_tag(connor._build_bam_tags(), 'X0')
        self.assertEqual('X0', tag._tag_name)
        self.assertEqual('Z', tag._tag_type)
        self.assertRegexpMatches(tag._description, 'filter')
        connor_align = MicroMock(filter='foo')
        self.assertEquals('foo', tag._get_value(None, connor_align))

    def test_build_bam_tags_x1_unique_identifier(self):
        tag = ConnorTest.get_tag(connor._build_bam_tags(), 'X1')
        self.assertEqual('X1', tag._tag_name)
        self.assertEqual('i', tag._tag_type)
        self.assertRegexpMatches(tag._description, 'unique identifier')
        family = MicroMock(umi_sequence=42)
        self.assertEquals(42, tag._get_value(family, None))

    def test_build_bam_tags_x2_umi_barcodes(self):
        tag = ConnorTest.get_tag(connor._build_bam_tags(), 'X2')
        self.assertEqual('X2', tag._tag_name)
        self.assertEqual('Z', tag._tag_type)
        self.assertRegexpMatches(tag._description, 'UMT barcodes')
        family = MicroMock(umi=('AAA','CCC'))
        self.assertEquals("AAA~CCC", tag._get_value(family, None))

    def test_build_bam_tags_x3_family_size(self):
        tag = ConnorTest.get_tag(connor._build_bam_tags(), 'X3')
        self.assertEqual('X3', tag._tag_name)
        self.assertEqual('i', tag._tag_type)
        self.assertRegexpMatches(tag._description, 'family size')
        family = MicroMock(alignments=[1]*42)
        self.assertEquals(42, tag._get_value(family, None))

    def test_build_bam_tags_x4_filter(self):
        tag = ConnorTest.get_tag(connor._build_bam_tags(), 'X4')
        self.assertEqual('X4', tag._tag_name)
        self.assertEqual('i', tag._tag_type)
        self.assertRegexpMatches(tag._description,
                                 'template for the consensus alignment')
        nontemplate_connor_align = MicroMock(query_name='foo')
        template_connor_align = MicroMock(query_name='bar')
        consensus = MicroMock(left_alignment=MicroMock(query_name='bar'))
        family = MicroMock(consensus=consensus)
        self.assertEquals(None,
                          tag._get_value(None, template_connor_align))
        self.assertEquals(None,
                          tag._get_value(family, nontemplate_connor_align))
        self.assertEquals(1, tag._get_value(family, template_connor_align))

    def test_build_annotated_aligns_writer(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>
'''.replace("|", "\t")

        with TempDirectory() as tmp_dir:
            input_bam = samtools_test.create_bam(tmp_dir.path,
                                                 'input.sam',
                                                 sam_contents)
            annotated_output_bam = os.path.join(tmp_dir.path, 'annotated.bam')
            tags = []
            actual_writer = connor._build_writer(input_bam,
                                                 annotated_output_bam,
                                                 tags)
            actual_writer.close()

            actual_output = samtools.alignment_file(annotated_output_bam, 'rb',)
            expected_header = {'HD': {'GO': 'none',
                                      'SO': 'coordinate',
                                      'VN': '1.4'},
                               'SQ': [{'SN': 'chr10', 'LN': 135534747}]}
            self.assertEqual(expected_header, actual_output.header)

    def test_build_annotated_aligns_writer_nullIfNotSpecified(self):
        actual_writer = connor._build_writer(input_bam='foo',
                                             output_bam='',
                                             tags=[])
        self.assertEqual(samtools.AlignWriter.NULL, actual_writer)

    def test_build_coordinate_read_name_manifest(self):
        Align = namedtuple('Align', 'name key')
        align1 = Align(name='align1', key=3)
        align2 = Align(name='align2', key=3)
        align3 = Align(name='align3', key=4)

        actual_dict = connor._build_coordinate_read_name_manifest([align1,
                                                           align2,
                                                           align3])

        expected_dict = {3: set(['align1', 'align2']),
                         4: set(['align3'])}
        self.assertEquals(expected_dict, actual_dict)

    def test_build_coordinate_families_oneFamily(self):
        align_A0 = MockAlignSegment("alignA", 'chr1', 10, 50, reference_end=16)
        align_A1 = MockAlignSegment("alignA", 'chr1', 50, 10, reference_end=56)
        align_B0 = MockAlignSegment("alignB", 'chr1', 10, 50, reference_end=16)
        align_B1 = MockAlignSegment("alignB", 'chr1', 50, 10, reference_end=56)
        alignments = [align_A0, align_B0, align_A1, align_B1]
        coord_read_name_manifest = {('chr1', 10, 56): set(['alignA', 'alignB'])}

        actual_families = [family for family in connor._build_coordinate_families(alignments, coord_read_name_manifest)]

        pair_A = connor._PairedAlignment(align_A0, align_A1)
        pair_B = connor._PairedAlignment(align_B0, align_B1)
        expected_families = [set([pair_A, pair_B])]
        self.assertEquals(expected_families, actual_families)

    def test_build_coordinate_families_threeFamilies(self):
        align_A0 = MockAlignSegment("alignA", 'chr1', 10, 70, reference_end=16)
        align_A1 = MockAlignSegment("alignA", 'chr1', 70, 10, reference_end=76)
        align_B0 = MockAlignSegment("alignB", 'chr1', 10, 70, reference_end=16)
        align_B1 = MockAlignSegment("alignB", 'chr1', 70, 10, reference_end=76)
        align_C0 = MockAlignSegment("alignC", 'chr1', 20, 80, reference_end=26)
        align_C1 = MockAlignSegment("alignC", 'chr1', 80, 20, reference_end=86)
        align_D0 = MockAlignSegment("alignD", 'chr1', 30, 90, reference_end=36)
        align_D1 = MockAlignSegment("alignD", 'chr1', 90, 30, reference_end=96)
        alignments = [align_A0, align_B0, align_C0, align_A1, align_B1,
                      align_D0, align_D1, align_C1]
        coord_read_name_manifest = {('chr1', 10, 76): set(['alignA', 'alignB']),
                                ('chr1', 20, 86): set(['alignC']),
                                ('chr1', 30, 96): set(['alignD'])}

        actual_families = [family for family in connor._build_coordinate_families(alignments,
                                                                                  coord_read_name_manifest)]

        pair_A = connor._PairedAlignment(align_A0, align_A1)
        pair_B = connor._PairedAlignment(align_B0, align_B1)
        pair_C = connor._PairedAlignment(align_C0, align_C1)
        pair_D = connor._PairedAlignment(align_D0, align_D1)

        expected_families = [set([pair_A, pair_B]),
                             set([pair_D]),
                             set([pair_C])]
        self.assertEquals(expected_families, actual_families)

    def test_build_tag_families_exact_left_or_right(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'NNNCCC')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'GGGNNN', 'NNNTTT')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'AAANNN', 'NNNCCC')
        pair4 = align_pair('alignD', 'chr1', 100, 200, 'TTTNNN', 'NNNCCC')
        pair5 = align_pair('alignE', 'chr1', 100, 200, 'CCCNNN', 'NNNTTT')

        paired_aligns = [pair1, pair2, pair3, pair4, pair5]
        tag1 = ('AAA', 'CCC')
        tag2 = ('GGG', 'TTT')
        ranked_tags = [tag1, tag2]

        actual_tag_fam_list = connor._build_tag_families(paired_aligns,
                                                         ranked_tags,
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)
        actual_tag_fam = {fam.umi:fam for fam in actual_tag_fam_list}

        self.assertEquals(2, len(actual_tag_fam))
        self.assertEquals(set([pair1, pair3, pair4]),
                          set(actual_tag_fam[tag1].alignments))
        self.assertEquals(set([pair2, pair5]),
                          set(actual_tag_fam[tag2].alignments))

    def test_build_tag_families_mostPopularTagIsCanonical(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'NNNCCC')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'AAANNN', 'NNNGGG')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'AAANNN', 'NNNGGG')
        pair4 = align_pair('alignD', 'chr1', 100, 200, 'TTTNNN', 'NNNGGG')
        pair5 = align_pair('alignE', 'chr1', 100, 200, 'TTTNNN', 'NNNCCC')
        input_pairs = [pair1, pair2, pair3, pair4, pair5]
        ranked_tags = [('AAA','GGG'),
                       ('AAA','CCC'),
                       ('TTT', 'CCC'),
                       ('TTT', 'GGG')]

        actual_tag_fam_list = connor._build_tag_families(input_pairs,
                                                         ranked_tags,
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)
        actual_tag_fam = {fam.umi:fam for fam in actual_tag_fam_list}

        self.assertEquals(set([pair1, pair2, pair3, pair4]),
                          set(actual_tag_fam[('AAA','GGG')].alignments))
        self.assertEquals(set([pair5]),
                          set(actual_tag_fam[('AAA','CCC')].alignments))

    def test_build_tag_families_exact_LR_included(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'NNNCCC')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(1, len(actual_tag_fam_list))
        self.assertEquals([pair1], actual_tag_fam_list[0].alignments)

    def test_build_tag_families_exact_L_included(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'NNNNNN')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(1, len(actual_tag_fam_list))
        self.assertEquals([pair1], actual_tag_fam_list[0].alignments)

    def test_build_tag_families_exact_R_included(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'NNNNNN', 'NNNCCC')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(1, len(actual_tag_fam_list))
        self.assertEquals([pair1], actual_tag_fam_list[0].alignments)

    def test_build_tag_families_hamming_L_included(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAGNNN', 'NNNNNN')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(1, len(actual_tag_fam_list))
        self.assertEquals([pair1], actual_tag_fam_list[0].alignments)

    def test_build_tag_families_hamming_R_included(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'NNNNNN', 'NNNCCG')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(1, len(actual_tag_fam_list))
        self.assertEquals([pair1], actual_tag_fam_list[0].alignments)

    def test_build_tag_families_hamming_LR_included(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAGNNN', 'NNNCCG')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(1, len(actual_tag_fam_list))
        self.assertEquals([pair1], actual_tag_fam_list[0].alignments)

    def test_build_tag_families_no_match_excluded(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'NNNNNN', 'NNNNNN')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(0, len(actual_tag_fam_list))

    def test_build_tag_families_hamming_dist_excluded(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AGGNNN', 'NNNCGG')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(0, len(actual_tag_fam_list))

    def test_build_tag_families_inexact_match_fuzzy_counted(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAGNNN', 'NNNNNN')

        families = connor._build_tag_families([pair1],
                                              [('AAA', 'CCC')],
                                              hamming_threshold=1,
                                              consensus_threshold=0.6)

        self.assertEquals(1, next(iter(families)).inexact_match_count)

    def test_build_tag_families_inexact_match_halfmatched_counted(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'NNNNNN')

        families = connor._build_tag_families([pair1],
                                              [('AAA', 'CCC')],
                                              hamming_threshold=1,
                                              consensus_threshold=0.6)

        self.assertEquals(1, next(iter(families)).inexact_match_count)

    def test_build_tag_families_inexact_match_not_counted(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'NNNCCC')

        families = connor._build_tag_families([pair1],
                                              [('AAA', 'CCC')],
                                              hamming_threshold=1,
                                              consensus_threshold=0.6)

        self.assertEquals(0, next(iter(families)).inexact_match_count)


    def test_hamming_distance_trivial(self):
        self.assertEqual(1, connor._hamming_dist("ABC", "ABD"))

    def test_hamming_distance_max(self):
        self.assertEqual(3, connor._hamming_dist("ABC", "XYZ"))

    def test_hamming_distance_0(self):
        self.assertEqual(0, connor._hamming_dist("ABC", "ABC"))

    def test_hamming_distance_unequal_lengths(self):
        self.assertRaises(AssertionError,
                          connor._hamming_dist,
                          "ABC",
                          "AB")

    def test_rank_tags_sortsByPopularity(self):
        pair0 = align_pair("align0", 'chr1', 100, 200, "TTTNNN", "NNNGGG")
        pair1 = align_pair("align1", 'chr1', 100, 200, "AAANNN", "NNNCCC")
        pair2 = align_pair("align2", 'chr1', 100, 200, "AAANNN", "NNNGGG")
        pair3 = align_pair("align3", 'chr1', 100, 200, "AAANNN", "NNNGGG")
        pair4 = align_pair("align4", 'chr1', 100, 200, "AAANNN", "NNNCCC")
        pair5 = align_pair("align5", 'chr1', 100, 200, "AAANNN", "NNNGGG")
        input_aligns = [pair0, pair1, pair2, pair3, pair4, pair5]

        actual_tags = connor._rank_tags(input_aligns)

        expected_tags = [('AAA', 'GGG'), ('AAA', 'CCC'), ('TTT', 'GGG')]
        self.assertEquals(expected_tags, actual_tags)

    def test_parse_command_line_args(self):
        namespace = connor._parse_command_line_args(["input.bam",
                                                     "output.bam"])
        self.assertEquals("input.bam", namespace.input_bam)
        self.assertEquals("output.bam", namespace.output_bam)

    def test_parse_command_line_args_throwsConnorUsageError(self):
        self.assertRaises(utils.UsageError,
                          connor._parse_command_line_args,
                          ["input"])
        self.assertRaises(utils.UsageError,
                          connor._parse_command_line_args,
                          ["input",
                           "output",
                           "something else"])


    def test_rank_tags_breaksTiesByTag(self):
        pair0 = align_pair("align0", 'chr1', 100, 200, "TTTNNN", "NNNGGG")
        pair1 = align_pair("align1", 'chr1', 100, 200, "AAANNN", "NNNCCC")
        pair2 = align_pair("align2", 'chr1', 100, 200, "AAANNN", "NNNGGG")
        input_aligns = [pair0, pair1, pair2]

        actual_tags = connor._rank_tags(input_aligns)

        expected_tags = [('AAA', 'CCC'), ('AAA', 'GGG'), ('TTT', 'GGG')]
        self.assertEquals(expected_tags, actual_tags)


class TestLightweightAlignment(BaseConnorTestCase):
    def test_lightweight_alignment_forwardRead(self):
        alignedSegment = align_seg("align1", 'chr1', 10, 100)

        actual_lwa = connor._LightweightAlignment(alignedSegment)

        self.assertEquals("align1", actual_lwa.name)
        self.assertEquals(('chr1', 10, 100), actual_lwa.key)

    def test_lightweight_alignment_reverseRead(self):
        alignedSegment = align_seg("align1", 'chr1', 100, 10)

        actual_lwa = connor._LightweightAlignment(alignedSegment)

        self.assertEquals("align1", actual_lwa.name)
        self.assertEquals(('chr1', 10, 100), actual_lwa.key)

    def test_lightweight_alignment_weirdRead(self):
        alignedSegment = align_seg("align1", 'chr1', 100, 100)

        actual_lwa = connor._LightweightAlignment(alignedSegment)

        self.assertEquals("align1", actual_lwa.name)
        self.assertEquals(('chr1', 100, 100), actual_lwa.key)

class ConnorIntegrationTestCase(BaseConnorTestCase):
    def test_deduplicate_alignnemnts(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>
readNameA2|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>
readNameB1|99|chr10|200|20|5M|=|400|200|CCCCC|>>>>>
readNameA1|147|chr10|300|20|5M|=|100|100|AAAAA|>>>>>
readNameA2|147|chr10|300|20|5M|=|100|100|AAAAA|>>>>>
readNameB1|147|chr10|400|20|5M|=|200|100|CCCCC|>>>>>
'''.replace("|", "\t")

        with TempDirectory() as tmp_dir:
            input_bam = samtools_test.create_bam(tmp_dir.path,
                                                 "input.sam",
                                                 sam_contents)
            output_bam = os.path.join(tmp_dir.path, "output.bam")
            consensus_writer = connor._build_writer(input_bam, output_bam, [])
            annotated_writer = samtools_test.MockAlignWriter()
            args = Namespace(input_bam=input_bam,
                             consensus_freq_threshold=0.6,
                             min_family_size_threshold=0,
                             umi_distance_threshold=1)
            connor._dedup_alignments(args,
                                     consensus_writer,
                                     annotated_writer,
                                     self.mock_logger)
            consensus_writer.close()
            alignments = samtools.alignment_file(output_bam, "rb").fetch()

            aligns = [(a.query_name, a.reference_start + 1) for a in alignments]
            self.assertEquals(4, len(aligns))
            self.assertEquals([("readNameA1", 100),
                               ("readNameB1", 200),
                               ("readNameA1", 300),
                               ("readNameB1", 400)],
                              aligns)

    def test_dedup_logging(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>
readNameA2|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>
readNameB1|99|chr10|200|20|5M|=|400|200|CCCCC|>>>>>
readNameA1|147|chr10|300|20|5M|=|100|100|AAAAA|>>>>>
readNameA2|147|chr10|300|20|5M|=|100|100|AAAAA|>>>>>
readNameB1|147|chr10|400|20|5M|=|200|100|CCCCC|>>>>>
'''.replace("|", "\t")

        with TempDirectory() as tmp_dir:
            input_bam = samtools_test.create_bam(tmp_dir.path,
                                                 'input.sam',
                                                 sam_contents)
            output_bam = os.path.join(tmp_dir.path, 'output.bam')
            args = Namespace(input_bam=input_bam,
                             output_bam=output_bam,
                             consensus_freq_threshold=0.6,
                             min_family_size_threshold=0,
                             umi_distance_threshold=1,
                             annotated_output_bam=None)
            connor._dedup_alignments(args,
                                     samtools_test.MockAlignWriter(),
                                     samtools_test.MockAlignWriter(),
                                     self.mock_logger)

            log_calls = self.mock_logger._log_calls['INFO']
            self.assertRegexpMatches(log_calls[0],
                                     'reading input bam')
            self.assertRegexpMatches(log_calls[1],
                                    'families were excluded')


    def test_deduplicate_alignments_distinctPairStartsAreNotCombined(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>
readNameB1|99|chr10|100|20|5M|=|500|200|AAAAA|>>>>>
readNameA1|147|chr10|300|20|5M|=|100|200|AAAAA|>>>>>
readNameB1|147|chr10|500|20|5M|=|100|200|AAAAA|>>>>>
'''.replace("|", "\t")

        with TempDirectory() as tmp_dir:
            input_bam = samtools_test.create_bam(tmp_dir.path,
                                                 "input.sam",
                                                 sam_contents)
            output_bam = os.path.join(tmp_dir.path, "output.bam")
            args = Namespace(input_bam=input_bam,
                             consensus_freq_threshold=0.6,
                             min_family_size_threshold=0,
                             umi_distance_threshold=1)
            consensus_writer = connor._build_writer(input_bam,
                                                    output_bam,
                                                    tags=[])
            annotated_writer = samtools.AlignWriter.NULL

            connor._dedup_alignments(args,
                                     consensus_writer,
                                     annotated_writer,
                                     self.mock_logger)
            consensus_writer.close()
            alignments = samtools_test.pysam_alignments_from_bam(output_bam)

            aligns = [(a.query_name, a.reference_start + 1) for a in alignments]
            self.assertEquals(4, len(aligns))
            self.assertEquals([("readNameA1", 100),
                               ("readNameB1", 100),
                               ("readNameA1", 300),
                               ("readNameB1", 500)],
                              aligns)


    def test_main_logging(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>
readNameA2|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>
readNameB1|99|chr10|200|20|5M|=|400|200|CCCCC|>>>>>
readNameA1|147|chr10|300|20|5M|=|100|100|AAAAA|>>>>>
readNameA2|147|chr10|300|20|5M|=|100|100|AAAAA|>>>>>
readNameB1|147|chr10|400|20|5M|=|200|100|CCCCC|>>>>>
'''.replace("|", "\t")

        with TempDirectory() as tmp_dir:
            input_bam = samtools_test.create_bam(tmp_dir.path,
                                                 'input.sam',
                                                 sam_contents)
            output_bam = os.path.join(tmp_dir.path, 'output.bam')
            output_log = os.path.join(tmp_dir.path, 'output.log')
            old_dedup_alignments = connor._dedup_alignments
            ##pylint: disable=unused-argument
            def angry_dedup(args, consensus_writer, annotated_writer, log):
                log.warning("possible problem")
            old_stderr = sys.stderr
            console_stream = StringIO()
            try:
                sys.stderr = console_stream
                connor._dedup_alignments = angry_dedup
                connor.main(["program_name",
                             input_bam,
                             output_bam,
                             "--min_family_size_threshold=0",
                             "--log_file=" + output_log])
            finally:
                connor._dedup_alignments = old_dedup_alignments
                sys.stderr = old_stderr
            log_lines = console_stream.getvalue().strip().split('\n')

        self.assertRegexpMatches(log_lines[0],
                                 r'connor begins \(v.*\)')
        self.assertRegexpMatches(log_lines[1],
                                 (r'logging to \[' + output_log + r'\]'))
        self.assertRegexpMatches(log_lines[2],
                                 r'possible problem')
        self.assertRegexpMatches(log_lines[3],
                                 (r'connor complete \(.*seconds.*memory\). '
                                 r'\*\*See warnings above\*\*'))
        self.assertEqual(4, len(log_lines))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    import unittest
    unittest.main()
