#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments
#pylint: disable=deprecated-method

from __future__ import print_function, absolute_import, division
from argparse import Namespace
from collections import namedtuple
import os
import sys
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from testfixtures.tempdirectory import TempDirectory
import connor.connor as connor
from connor import samtools
import connor.utils as utils
from connor.connor import _build_lightweight_pairs
from connor.samtools import ConnorAlign
from connor.connor import _PairedAlignment
from test.samtools_test import MockAlignWriter
from test.samtools_test import mock_align
import test.samtools_test as samtools_test
from test.utils_test import BaseConnorTestCase
from test.utils_test import MicroMock


def _mock_connor_align(query_name,
                       reference_name,
                       reference_start,
                       next_reference_start,
                       reference_end):
    mock_pysam = MicroMock(query_name=query_name,
                           reference_name=reference_name,
                           reference_start=reference_start,
                           next_reference_start=next_reference_start,
                           reference_end=reference_end,
                           query_sequence='AAAGGG')
    return ConnorAlign(mock_pysam)




def _mock_tag_family(align_pairs=None,
                    distinct_cigar_count=1,
                    inexact_match_count=0,
                    minority_cigar_percentage=0,
                    consensus=None,
                    filter_value=None,
                    included_pair_count=5):
    if align_pairs is None:
        align_pairs = [1,2,3,4]
    return MicroMock(align_pairs=align_pairs,
                     distinct_cigar_count=distinct_cigar_count,
                     inexact_match_count=inexact_match_count,
                     minority_cigar_percentage=minority_cigar_percentage,
                     consensus=consensus,
                     filter_value=filter_value,
                     included_pair_count=included_pair_count)


# #TODO: cgates: replace this with samtools_test.mock_align
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
        self.filter_value = None
        self.mapping_quality = 20

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

        self.assertIs(left_align, actual_paired_alignment.left)
        self.assertIs(right_align, actual_paired_alignment.right)
        self.assertEquals(("AAATTT","CCCGGG"), actual_paired_alignment.umt)

    def test_init_valueErrorOnInconsistentQueryNames(self):
        left = mock_align(query_name="alignA")
        right = mock_align(query_name="alignB")
        self.assertRaisesRegexp(ValueError,
                                r'Inconsistent query names \(alignA != alignB\)',
                                connor._PairedAlignment,
                                left,
                                right,
                                tag_length=1)

    def test_filter_value(self):
        left = ConnorAlign(mock_align(), filter_value=None)
        right = ConnorAlign(mock_align(), filter_value=None)
        paired_alignment = connor._PairedAlignment(left, right, tag_length=1)
        self.assertEqual(None, paired_alignment.filter_value)

        left = ConnorAlign(mock_align(), filter_value='')
        right = ConnorAlign(mock_align(), filter_value='')
        paired_alignment = connor._PairedAlignment(left, right, tag_length=1)
        self.assertEqual(None, paired_alignment.filter_value)

        left = ConnorAlign(mock_align(), filter_value='foo')
        right = ConnorAlign(mock_align(), filter_value=None)
        paired_alignment = connor._PairedAlignment(left, right, tag_length=1)
        self.assertEqual(('foo', None), paired_alignment.filter_value)

        left = ConnorAlign(mock_align(), filter_value=None)
        right = ConnorAlign(mock_align(), filter_value='bar')
        paired_alignment = connor._PairedAlignment(left, right, tag_length=1)
        self.assertEqual((None, 'bar'), paired_alignment.filter_value)

    def test_query_name(self):
        left = mock_align(query_name="alignA", reference_start=100)
        right = mock_align(query_name="alignA", reference_start=200)
        paired_alignment = connor._PairedAlignment(left, right, tag_length=1)
        self.assertEqual("alignA", paired_alignment.query_name)

    def test_eq(self):
        left = align_seg("alignA", 'chr1', 100, 200, "AAANNNNNNN")
        right = align_seg("alignA", 'chr1', 200, 100, "CCCNNNNNNN")
        other = align_seg("alignA", 'chr2', 100, 200, "AAANNNNNNN")

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

    def test_replace_umt(self):
        # pysam's represtation of the sequence is inconsistent across pysam
        # and python versions; this hack makes the values comparable
        def _byte_array_to_string(sequence):
            if isinstance(sequence, str):
                return sequence
            else:
                return str(sequence.decode("utf-8"))

        left_A = mock_align(query_sequence='AANN', query_qualities=[1,2,3,4])
        right_A = mock_align(query_sequence='NNCC', query_qualities=[5,6,7,8])
        paired_align = connor._PairedAlignment(left_A, right_A, tag_length=2)

        paired_align.replace_umt(('GG','TT'))

        self.assertEquals('GGNN',
                          _byte_array_to_string(paired_align.left.query_sequence))
        self.assertEquals('NNTT',
                          _byte_array_to_string(paired_align.right.query_sequence))
        self.assertEquals([1,2,3,4],
                          paired_align.left.query_qualities)
        self.assertEquals([5,6,7,8],
                          paired_align.right.query_qualities)

    def test_replace_umt_errorIfInconsistentUmtLength(self):
        left_A = mock_align(query_sequence='AANN', query_qualities=[1,2,3,4])
        right_A = mock_align(query_sequence='NNCC', query_qualities=[5,6,7,8])
        paired_align = connor._PairedAlignment(left_A, right_A, tag_length=2)

        self.assertRaisesRegexp(ValueError,
                                r'Each UMT must match tag_length \(2\)',
                                paired_align.replace_umt,
                                ('G','TT'))
        self.assertRaisesRegexp(ValueError,
                                r'Each UMT must match tag_length \(2\)',
                                paired_align.replace_umt,
                                ('GG','T'))
        self.assertRaisesRegexp(ValueError,
                                r'Each UMT must match tag_length \(2\)',
                                paired_align.replace_umt,
                                (None, None))
        self.assertRaisesRegexp(ValueError,
                                r'Each UMT must match tag_length \(2\)',
                                paired_align.replace_umt,
                                ('G',))


class TagFamiliyTest(BaseConnorTestCase):
    def test_init(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'GGGNNN', 'TTTNNN')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        align_pairs = [pair1, pair2, pair3]
        inexact_match_count = 42

        input_umt = "AAANNNCCCNNN"
        actual_tag_family = connor._TagFamily(input_umt,
                                             align_pairs,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(input_umt, actual_tag_family.umt)
        self.assertEquals(align_pairs, actual_tag_family.align_pairs)
        self.assertEquals(42, actual_tag_family.inexact_match_count)

    def test_init_setsFilterValue(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        align_pairs = [pair1, pair2, pair3]
        inexact_match_count = 42
        input_umt = "AAANNNCCCNNN"

        def family_filter(family):
            return ":".join([a.query_name for a in family.align_pairs]) +\
                 " is too odd"
        actual_tag_family = connor._TagFamily(input_umt,
                                              align_pairs,
                                              inexact_match_count,
                                              consensus_threshold=0.6,
                                              family_filter=family_filter)

        self.assertEquals('alignA:alignB:alignC is too odd',
                          actual_tag_family.filter_value)


    def test_consensus_rewrites_umt(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'GGGNNN', 'NNNTTT')
        input_umts = ("AAA", "CCC")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umts,
                                             [pair1],
                                             inexact_match_count,
                                             consensus_threshold=0.6)
        actual_consensus = actual_tag_family.consensus

        self.assertEquals(input_umts, actual_consensus.umt)
        self.assertEquals('AAANNN',
                          actual_consensus.left.query_sequence)
        self.assertEquals('NNNCCC',
                          actual_consensus.right.query_sequence)

    def test_consensus_sequence_trivial_noop(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)
        actual_consensus_seq = actual_tag_family.consensus

        self.assertEquals("nnnGGG",
                          actual_consensus_seq.left.query_sequence)
        self.assertEquals("TTTnnn",
                          actual_consensus_seq.right.query_sequence)

    def test_consensus_sequence_majority_wins(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'TCTnnn')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'TCTnnn')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")

        actual_tag_family = connor._TagFamily(input_umts,
                                              alignments,
                                              inexact_match_count=1,
                                              consensus_threshold=0.5)
        consensus_pair = actual_tag_family.consensus

        self.assertEquals("nnnGTG",
                          consensus_pair.left.query_sequence)
        self.assertEquals("TCTnnn",
                          consensus_pair.right.query_sequence)

    def test_consensus_sequence_below_threshold_Ns(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'TCTnnn')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'TCTnnn')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        threshold = 0.8
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             threshold)
        consensus_pair = actual_tag_family.consensus

        self.assertEquals("nnnGNG",
                          consensus_pair.left.query_sequence)
        self.assertEquals("TNTnnn",
                          consensus_pair.right.query_sequence)

    def connor_align(self, query_name, query_sequence, mapping_quality):
        return ConnorAlign(mock_align(query_name=query_name,
                                      query_sequence=query_sequence,
                                      mapping_quality=mapping_quality))

    def test_consensus_qualities_maxMappingQualityScores(self):
        alignAL = self.connor_align('alignA', 'nGT', 30)
        alignAR = self.connor_align('alignA', 'nCT', 25)
        pairA = connor._PairedAlignment(alignAL, alignAR, tag_length=1)

        alignBL = self.connor_align('alignB', 'nGT', 20)
        alignBR = self.connor_align('alignB', 'nCT', 15)
        pairB = connor._PairedAlignment(alignBL, alignBR, tag_length=1)

        alignCL = self.connor_align('alignC', 'nGT', 10)
        alignCR = self.connor_align('alignC', 'nCT', 5)
        pairC = connor._PairedAlignment(alignCL, alignCR, tag_length=1)

        alignments = [pairA, pairB, pairC]
        input_umts = ("n", "n")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umts,
                                               alignments,
                                               inexact_match_count,
                                               consensus_threshold=0.6)
        paired_consensus = actual_tag_family.consensus

        self.assertEquals(30,
                          paired_consensus.left.mapping_quality)
        self.assertEquals(25,
                          paired_consensus.right.mapping_quality)


    def test_select_template_alignment_pair_picksMaxQualityScores(self):
        alignAL = self.connor_align('alignA', 'nGT', 20)
        alignAR = self.connor_align('alignA', 'nCT', 15)
        pairA = connor._PairedAlignment(alignAL, alignAR, tag_length=1)

        alignBL = self.connor_align('alignB', "nGT", 30)
        alignBR = self.connor_align('alignB', "nCT", 25)
        pairB = connor._PairedAlignment(alignBL, alignBR, tag_length=1)

        alignCL = self.connor_align('alignC', "nGT", 10)
        alignCR = self.connor_align('alignC', "nCT", 5)
        pairC = connor._PairedAlignment(alignCL, alignCR, tag_length=1)

        alignment_pairs = [pairA, pairB, pairC]

        actual_template = connor._TagFamily._select_template_alignment_pair(alignment_pairs)

        self.assertEquals(pairB, actual_template)

    def test_select_template_alignment_pair_breaksTiesByQueryName(self):
        alignAL = self.connor_align('alignA', "nGT", 20)
        alignAR = self.connor_align('alignA', "nCT", 15)
        pairA = connor._PairedAlignment(alignAL, alignAR, tag_length=1)

        alignBL = self.connor_align('alignB', "nGT", 20)
        alignBR = self.connor_align('alignB', "nCT", 15)
        pairB = connor._PairedAlignment(alignBL, alignBR, tag_length=1)
        alignment_pairs = [pairA, pairB]

        actual_template = connor._TagFamily._select_template_alignment_pair(alignment_pairs)

        self.assertEquals(pairA, actual_template)


    def test_consensus_uniform_cigars_admitted(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left.cigarstring = "3S3M"
        pair1.right.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left.cigarstring = "3S3M"
        pair2.right.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left.cigarstring = "3S3M"
        pair3.right.cigarstring = "3S3M"
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        included_names = [a.left.query_name for a in actual_tag_family.align_pairs]
        self.assertEquals(['alignA','alignB','alignC'], included_names)

    def test_consensus_minority_cigars_excluded(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left.cigarstring = "3S3M"
        pair1.right.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left.cigarstring = "3S3M"
        pair2.right.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left.cigarstring = "3S3I"
        pair3.right.cigarstring = "3S3M"
        align_pairs = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umts,
                                              align_pairs,
                                              inexact_match_count,
                                              consensus_threshold=0.6)

        pair_name_filter = {}
        for pair in actual_tag_family.align_pairs:
            query_name = pair.left.query_name
            pair_name_filter[query_name] = (pair.left.filter_value,
                                            pair.right.filter_value)
        self.assertEquals((None, None), pair_name_filter['alignA'])
        self.assertEquals((None, None), pair_name_filter['alignB'])
        self.assertEquals(('minority CIGAR', 'minority CIGAR'),
                          pair_name_filter['alignC'])


    def test_init_distinctCigarCountTwo(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left.cigarstring = "3S3M"
        pair1.right.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left.cigarstring = "3S3M"
        pair2.right.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left.cigarstring = "3S3I"
        pair3.right.cigarstring = "3S3I"
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(2, actual_tag_family.distinct_cigar_count)

    def test_init_distinctCigarCountOneForConsisentCigar(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left.cigarstring = "3S3M"
        pair1.right.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left.cigarstring = "3S3M"
        pair2.right.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left.cigarstring = "3S3M"
        pair3.right.cigarstring = "3S3M"
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(1, actual_tag_family.distinct_cigar_count)

    def test_distinct_cigar_count_single_ended(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left.cigarstring = "3S3M"
        pair1.right.cigarstring = "3S3M"
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left.cigarstring = "3S3M"
        pair2.right.cigarstring = "3S3M"
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left.cigarstring = "3S3I"
        pair3.right.cigarstring = "3S3M"
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(2, actual_tag_family.distinct_cigar_count)

#     def test_init_noMinorityIfConsistentCigar(self):
#         pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
#         pair1.left.cigarstring = "3S3M"
#         pair1.right.cigarstring = "3S3M"
#         pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
#         pair2.left.cigarstring = "3S3M"
#         pair2.right.cigarstring = "3S3M"
#         pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
#         pair3.left.cigarstring = "3S3M"
#         pair3.right.cigarstring = "3S3M"
#         alignments = [pair1, pair2, pair3]
#         input_umts = ("nnn", "nnn")
#         inexact_match_count = 0
# 
#         actual_tag_family = connor._TagFamily(input_umts,
#                                              alignments,
#                                              inexact_match_count,
#                                              consensus_threshold=0.6)
# 
#         self.assertEquals(0, actual_tag_family.minority_cigar_percentage)
# 
#     def test_init_setsMinorityPercentage(self):
#         pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
#         pair1.left.cigarstring = "3S3M"
#         pair1.right.cigarstring = "3S3M"
#         pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
#         pair2.left.cigarstring = "3S3M"
#         pair2.right.cigarstring = "3S3M"
#         pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
#         pair3.left.cigarstring = "3S3I"
#         pair3.right.cigarstring = "3S3I"
#         alignments = [pair1, pair2, pair3]
#         input_umts = ("nnn", "nnn")
#         inexact_match_count = 0
# 
#         actual_tag_family = connor._TagFamily(input_umts,
#                                              alignments,
#                                              inexact_match_count,
#                                              consensus_threshold=0.6)
# 
#         self.assertEquals(1/3, actual_tag_family.minority_cigar_percentage)

    def test_init_addsFilterForMinority(self):
        def pair(name, c1, c2):
            return _PairedAlignment(ConnorAlign(mock_align(query_name=name,
                                                           cigarstring=c1)),
                                    ConnorAlign(mock_align(query_name=name,
                                                           cigarstring=c2)))
        pairA = pair('alignA', '3S3M', '3S3M')
        pairB = pair('alignB', '3S3M', '3S3M')
        pairC = pair('alignC', '3S1I3M', '3S3M')
        align_pairs = [pairA, pairB, pairC]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umts,
                                             align_pairs,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        pair_name_filter = {}
        for pair in actual_tag_family.align_pairs:
            query_name = pair.left.query_name
            pair_name_filter[query_name] = (pair.left.filter_value,
                                            pair.right.filter_value)

        self.assertEqual(3, len(pair_name_filter))
        self.assertEqual((None, None), pair_name_filter['alignA'])
        self.assertEqual((None, None), pair_name_filter['alignB'])
        self.assertEqual(('minority CIGAR','minority CIGAR'),
                         pair_name_filter['alignC'])


    def test_included_pair_count(self):
        def pair(c1, c2):
            return _PairedAlignment(ConnorAlign(mock_align(cigarstring=c1)),
                                    ConnorAlign(mock_align(cigarstring=c2)))
        pairA = pair("3S3M", "3S3M")
        pairB = pair("3S3M", "3S3M")
        pairC = pair("3S1I3M", "3S3M")
        alignments = [pairA, pairB, pairC]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = connor._TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(2, actual_tag_family.included_pair_count)


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

        self.assertEquals(None, tag._get_value(None, None))

        family = MicroMock(filter_value=None)
        connor_align = MicroMock(filter_value=None)
        self.assertEquals(None, tag._get_value(family, connor_align))

        family = MicroMock(filter_value='foo')
        connor_align = MicroMock(filter_value='bar')
        self.assertEquals('foo', tag._get_value(family, None))
        self.assertEquals('bar', tag._get_value(None, connor_align))
        self.assertEquals('foo; bar', tag._get_value(family, connor_align))



    def test_build_bam_tags_x1_unique_identifier(self):
        tag = ConnorTest.get_tag(connor._build_bam_tags(), 'X1')
        self.assertEqual('X1', tag._tag_name)
        self.assertEqual('i', tag._tag_type)
        self.assertRegexpMatches(tag._description, 'unique identifier')
        family = MicroMock(umi_sequence=42)
        self.assertEquals(42, tag._get_value(family, None))
        self.assertEquals(None, tag._get_value(None, None))


    def test_build_bam_tags_x2_umt_barcodes(self):
        tag = ConnorTest.get_tag(connor._build_bam_tags(), 'X2')
        self.assertEqual('X2', tag._tag_name)
        self.assertEqual('Z', tag._tag_type)
        self.assertRegexpMatches(tag._description, 'UMT barcodes')
        family = MicroMock(umt=('AAA','CCC'))
        self.assertEquals("AAA~CCC", tag._get_value(family, None))
        self.assertEquals(None, tag._get_value(None, None))


    def test_build_bam_tags_x3_family_size(self):
        tag = ConnorTest.get_tag(connor._build_bam_tags(), 'X3')
        self.assertEqual('X3', tag._tag_name)
        self.assertEqual('i', tag._tag_type)
        self.assertRegexpMatches(tag._description, 'family size')
        family = MicroMock(included_pair_count=42)
        self.assertEquals(42, tag._get_value(family, None))
        self.assertEquals(None, tag._get_value(None, None))

    def test_build_bam_tags_x4_filter(self):
        tag = ConnorTest.get_tag(connor._build_bam_tags(), 'X4')
        self.assertEqual('X4', tag._tag_name)
        self.assertEqual('i', tag._tag_type)
        self.assertRegexpMatches(tag._description,
                                 'template for the consensus alignment')
        nontemplate_connor_align = MicroMock(query_name='foo')
        template_connor_align = MicroMock(query_name='bar')
        consensus_pair = MicroMock(left=MicroMock(query_name='bar'))
        family = MicroMock(consensus=consensus_pair)
        self.assertEquals(None,
                          tag._get_value(None, template_connor_align))
        self.assertEquals(None,
                          tag._get_value(family, nontemplate_connor_align))
        self.assertEquals(1, tag._get_value(family, template_connor_align))
        self.assertEquals(None, tag._get_value(None, None))


    def test_build_family_filter_whenFamilySizeOk(self):
        args = Namespace(min_family_size_threshold=2)
        family_filter = connor._build_family_filter(args)
        family = MicroMock(included_pair_count=2)
        self.assertEqual(None, family_filter(family))

    def test_build_family_filter_whenFamilySizeTooSmall(self):
        args = Namespace(min_family_size_threshold=3)
        family_filter = connor._build_family_filter(args)
        family = MicroMock(included_pair_count = 2)
        self.assertEqual('family too small (<3)', family_filter(family))

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
        align_A0 = _mock_connor_align("alignA", 'chr1', 10, 50, 16)
        align_A1 = _mock_connor_align("alignA", 'chr1', 50, 10, 56)
        align_B0 = _mock_connor_align("alignB", 'chr1', 10, 50, 16)
        align_B1 = _mock_connor_align("alignB", 'chr1', 50, 10, 56)
        alignments = [align_A0, align_B0, align_A1, align_B1]
        coord_read_name_manifest = {('chr1', 10, 56): set(['alignA', 'alignB'])}

        family_gen = connor._build_coordinate_families(alignments,
                                                       coord_read_name_manifest,
                                                       MockAlignWriter())
        actual_families = [family for family in family_gen]

        pair_A = connor._PairedAlignment(align_A0, align_A1)
        pair_B = connor._PairedAlignment(align_B0, align_B1)
        expected_families = [set([pair_A, pair_B])]
        self.assertEquals(expected_families, actual_families)

    def test_build_coordinate_families_writeRemainders(self):
        align_A0 = _mock_connor_align('alignA', 'chr1', 10, 50, 16)
        align_B0 = _mock_connor_align('alignB', 'chr1', 10, 50, 16)
        align_B1 = _mock_connor_align('alignB', 'chr1', 50, 10, 56)
        alignments = [align_A0, align_B0, align_B1]
        coord_read_name_manifest = {('chr1', 10, 56): set(['alignA', 'alignB'])}

        excluded_writer = MockAlignWriter()

        family_gen = connor._build_coordinate_families(alignments,
                                                       coord_read_name_manifest,
                                                       excluded_writer)
        for _dummy in family_gen:
            pass

        self.assertEqual(1, len(excluded_writer._write_calls))
        self.assertEqual((None, align_A0), excluded_writer._write_calls[0])

    def test_build_coordinate_families_writeRemaindersIsOrdered(self):
        align_C = _mock_connor_align('alignC', 'chr1', 10, 50, 16)
        align_A = _mock_connor_align('alignA', 'chr1', 10, 50, 16)
        align_B = _mock_connor_align('alignB', 'chr1', 50, 10, 56)
        alignments = [align_C, align_A, align_B]
        coord_read_name_manifest = {('chr1', 10, 56): set([])}

        excluded_writer = MockAlignWriter()

        family_gen = connor._build_coordinate_families(alignments,
                                                       coord_read_name_manifest,
                                                       excluded_writer)
        for _dummy in family_gen:
            pass

        excluded_names = [align.query_name for (_, align) in excluded_writer._write_calls]
        self.assertEqual(['alignA', 'alignB', 'alignC'],
                         excluded_names)


    def test_build_coordinate_families_threeFamilies(self):
        align_A0 = _mock_connor_align("alignA", 'chr1', 10, 70, 16)
        align_A1 = _mock_connor_align("alignA", 'chr1', 70, 10, 76)
        align_B0 = _mock_connor_align("alignB", 'chr1', 10, 70, 16)
        align_B1 = _mock_connor_align("alignB", 'chr1', 70, 10, 76)
        align_C0 = _mock_connor_align("alignC", 'chr1', 20, 80, 26)
        align_C1 = _mock_connor_align("alignC", 'chr1', 80, 20, 86)
        align_D0 = _mock_connor_align("alignD", 'chr1', 30, 90, 36)
        align_D1 = _mock_connor_align("alignD", 'chr1', 90, 30, 96)
        alignments = [align_A0, align_B0, align_C0, align_A1, align_B1,
                      align_D0, align_D1, align_C1]
        coord_read_name_manifest = {('chr1', 10, 76): set(['alignA', 'alignB']),
                                ('chr1', 20, 86): set(['alignC']),
                                ('chr1', 30, 96): set(['alignD'])}

        family_gen = connor._build_coordinate_families(alignments,
                                                       coord_read_name_manifest,
                                                       MockAlignWriter())
        actual_families = [family for family in family_gen]

        pair_A = connor._PairedAlignment(align_A0, align_A1)
        pair_B = connor._PairedAlignment(align_B0, align_B1)
        pair_C = connor._PairedAlignment(align_C0, align_C1)
        pair_D = connor._PairedAlignment(align_D0, align_D1)

        expected_families = [set([pair_A, pair_B]),
                             set([pair_D]),
                             set([pair_C])]
        self.assertEquals(expected_families, actual_families)

    def test_build_lightweight_alignment_pairs_simple(self):
        alignAL = mock_align(query_name='alignA', pos=100, query_sequence='AAA',
                             cigarstring='3M')
        alignAR = mock_align(query_name='alignA', pos=200, query_sequence='AAA',
                             cigarstring='3M')
        alignments = [alignAL, alignAR]

        lightweight_pairs = _build_lightweight_pairs(alignments)

        self.assertEqual(1, len(lightweight_pairs))
        self.assertEqual('alignA', lightweight_pairs[0].name)
        chrom = None
        left_pos_start = 100
        right_pos_end = 203
        self.assertEqual((chrom, left_pos_start, right_pos_end),
                         lightweight_pairs[0].key)

    def test_build_lightweight_alignment_pairs_singletonsDiscarded(self):
        alignAL = mock_align(query_name='alignA', pos=100)
        alignAR = mock_align(query_name='alignA', pos=200)
        alignBL = mock_align(query_name='alignB', pos=100)
        alignBR = mock_align(query_name='alignB', pos=200)
        alignCR = mock_align(query_name='alignC', pos=100)

        alignments = [alignAL, alignAR, alignBL, alignBR, alignCR]

        lightweight_pairs = _build_lightweight_pairs(alignments)

        self.assertEqual(2, len(lightweight_pairs))
        align_names = set([p.name for p in lightweight_pairs])
        self.assertTrue('alignA' in align_names)
        self.assertTrue('alignB' in align_names)

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
        actual_tag_fam = {fam.umt:fam for fam in actual_tag_fam_list}

        self.assertEquals(2, len(actual_tag_fam))
        self.assertEquals(set([pair1, pair3, pair4]),
                          set(actual_tag_fam[tag1].align_pairs))
        self.assertEquals(set([pair2, pair5]),
                          set(actual_tag_fam[tag2].align_pairs))

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
        actual_tag_fam = {fam.umt:fam for fam in actual_tag_fam_list}

        self.assertEquals(set([pair1, pair2, pair3, pair4]),
                          set(actual_tag_fam[('AAA','GGG')].align_pairs))
        self.assertEquals(set([pair5]),
                          set(actual_tag_fam[('AAA','CCC')].align_pairs))

    def test_build_tag_families_exact_LR_included(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'NNNCCC')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(1, len(actual_tag_fam_list))
        self.assertEquals(set([pair1]), actual_tag_fam_list[0].align_pairs)

    def test_build_tag_families_exact_L_included(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'NNNNNN')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(1, len(actual_tag_fam_list))
        self.assertEquals(set([pair1]), actual_tag_fam_list[0].align_pairs)

    def test_build_tag_families_exact_R_included(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'NNNNNN', 'NNNCCC')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(1, len(actual_tag_fam_list))
        self.assertEquals(set([pair1]), actual_tag_fam_list[0].align_pairs)

    def test_build_tag_families_hamming_L_included(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAGNNN', 'NNNNNN')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(1, len(actual_tag_fam_list))
        self.assertEquals(set([pair1]), actual_tag_fam_list[0].align_pairs)

    def test_build_tag_families_hamming_R_included(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'NNNNNN', 'NNNCCG')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(1, len(actual_tag_fam_list))
        self.assertEquals(set([pair1]), actual_tag_fam_list[0].align_pairs)

    def test_build_tag_families_hamming_LR_included(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAGNNN', 'NNNCCG')

        actual_tag_fam_list = connor._build_tag_families([pair1],
                                                         [('AAA', 'CCC')],
                                                         hamming_threshold=1,
                                                         consensus_threshold=0.6)

        self.assertEquals(1, len(actual_tag_fam_list))
        self.assertEquals(set([pair1]), actual_tag_fam_list[0].align_pairs)

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
        self.assertEqual(expected_tags, actual_tags)

    def test_parse_command_line_args(self):
        namespace = connor._parse_command_line_args(["command",
                                                     "input.bam",
                                                     "output.bam"])
        self.assertEqual("input.bam", namespace.input_bam)
        self.assertEqual("output.bam", namespace.output_bam)
        self.assertEqual(False, namespace.force)
        self.assertEqual(False, namespace.simplify_pg_header)
        self.assertEqual(False, namespace.verbose)
        self.assertEqual("output.bam.log", namespace.log_file)
        self.assertEqual(None, namespace.annotated_output_bam)
        self.assertEqual(connor.DEFAULT_CONSENSUS_FREQ_THRESHOLD,
                         namespace.consensus_freq_threshold)
        self.assertEqual(connor.DEFAULT_MIN_FAMILY_SIZE_THRESHOLD,
                         namespace.min_family_size_threshold)
        self.assertEqual(connor.DEFAULT_UMT_DISTANCE_THRESHOLD,
                         namespace.umt_distance_threshold)
        self.assertEquals(['command', 'input.bam', 'output.bam'],
                          namespace.original_command_line)
        self.assertEqual(11, len(vars(namespace)))

    def test_parse_command_line_args_throwsConnorUsageError(self):
        self.assertRaises(utils.UsageError,
                          connor._parse_command_line_args,
                          ["command", "input"])
        self.assertRaises(utils.UsageError,
                          connor._parse_command_line_args,
                          ["command",
                           "input",
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
        alignedSegment = _mock_connor_align("align1", 'chr1', 10, 100, 16)

        actual_lwa = connor._LightweightAlignment(alignedSegment)

        self.assertEquals("align1", actual_lwa.name)
        self.assertEquals(('chr1', 10, 100), actual_lwa.key)

    def test_lightweight_alignment_reverseRead(self):
        alignedSegment = _mock_connor_align("align1", 'chr1', 100, 10, 16)

        actual_lwa = connor._LightweightAlignment(alignedSegment)

        self.assertEquals("align1", actual_lwa.name)
        self.assertEquals(('chr1', 10, 100), actual_lwa.key)

    def test_lightweight_alignment_weirdRead(self):
        alignedSegment = _mock_connor_align("align1", 'chr1', 100, 100, 16)

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
            args = Namespace(simplify_pg_header=True,
                             original_command_line='foo')
            consensus_writer = samtools.build_writer(input_bam,
                                                     output_bam,
                                                     tags=[],
                                                     args=args)
            annotated_writer = samtools_test.MockAlignWriter()
            args = Namespace(input_bam=input_bam,
                             consensus_freq_threshold=0.6,
                             min_family_size_threshold=0,
                             umt_distance_threshold=1)
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
readNameC1|99|chr10|400|20|5M|=|200|100|CCCCC|>>>>>
readNameC1|147|chr10|400|20|5M|=|200|100|CCCCC|>>>>>
readNameC2|99|chr10|400|20|5M|=|200|100|CCCCC|>>>>>
readNameC2|147|chr10|400|20|5M|=|200|100|CCCCC|>>>>>
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
                             umt_distance_threshold=1,
                             annotated_output_bam=None)
            connor._dedup_alignments(args,
                                     samtools_test.MockAlignWriter(),
                                     samtools_test.MockAlignWriter(),
                                     self.mock_logger)

            log_iter = iter(self.mock_logger._log_calls['INFO'])
            self.assertRegexpMatches(next(log_iter),
                                     'reading input bam')
            self.assertRegexpMatches(next(log_iter),
                                     '0% .* alignments processed')
            self.assertRegexpMatches(next(log_iter),
                                     '20% .* alignments processed')
            self.assertRegexpMatches(next(log_iter),
                                     '30% .* alignments processed')
            self.assertRegexpMatches(next(log_iter),
                                     '40% .* alignments processed')
            self.assertRegexpMatches(next(log_iter),
                                     '50% .* alignments processed')
            self.assertRegexpMatches(next(log_iter),
                                     '60% .* alignments processed')
            self.assertRegexpMatches(next(log_iter),
                                     '70% .* alignments processed')
            self.assertRegexpMatches(next(log_iter),
                                     '80% .* alignments processed')
            self.assertRegexpMatches(next(log_iter),
                                     '90% .* alignments processed')
            self.assertRegexpMatches(next(log_iter),
                                     '100% .* alignments processed')


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
                             umt_distance_threshold=1,
                             simplify_pg_header=False,
                             original_command_line='foo')
            consensus_writer = samtools.build_writer(input_bam,
                                                     output_bam,
                                                     [],
                                                     args)
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
            #pylint: disable=unused-argument
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
        self.assertRegexpMatches(log_lines[3], 'sorting/indexing')
        self.assertRegexpMatches(log_lines[4],
                                 (r'connor complete \(.*seconds.*memory\). '
                                 r'\*\*See warnings above\*\*'))
        self.assertEqual(5, len(log_lines))


if __name__ == "__main__":
    import unittest
    unittest.main()
