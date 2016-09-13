#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments
#pylint: disable=deprecated-method

from __future__ import print_function, absolute_import, division
from argparse import Namespace
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
from connor.connor import _CoordinateFamilyHolder
from connor.samtools import ConnorAlign
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
    return samtools.PairedAlignment(alignL, alignR, tag_length)


class PairedAlignmentTest(BaseConnorTestCase):
    def test_init(self):
        left_align = MockAlignSegment("alignA", 'chr1', 10, 20,
                                      "AAATTT" "GGGG", reference_end=14)
        right_align = MockAlignSegment("alignA", 'chr1', 20, 10,
                                       "TTTT" "CCCGGG" ,reference_end=24)
        tag_length = 6
        actual_paired_alignment = samtools.PairedAlignment(left_align,
                                                         right_align,
                                                         tag_length)

        self.assertIs(left_align, actual_paired_alignment.left)
        self.assertIs(right_align, actual_paired_alignment.right)
        self.assertEquals(("AAATTT","CCCGGG"), actual_paired_alignment.umt)

    def test_init_valueErrorOnInconsistentQueryNames(self):
        left = mock_align(query_name="alignA")
        right = mock_align(query_name="alignB")
        self.assertRaisesRegexp(ValueError,
                                (r'Inconsistent query names '
                                 r'\(alignA != alignB\)'),
                                samtools.PairedAlignment,
                                left,
                                right,
                                tag_length=1)

    def test_cigars(self):
        left = MicroMock(query_name='A',
                         cigarstring='1S2M4S',
                         query_sequence='AAAAAA')
        right = MicroMock(query_name='A',
                         cigarstring='16S32M64S',
                          query_sequence='AAAAAA')
        paired_alignment = samtools.PairedAlignment(left, right, tag_length=1)
        self.assertEqual(('1S2M4S', '16S32M64S'), paired_alignment.cigars)

    def test_positions(self):
        left = MicroMock(query_name='A',
                         reference_start=100,
                         reference_end=150,
                         query_sequence='AAAAAA')
        right = MicroMock(query_name='A',
                          reference_start=200,
                          reference_end=250,
                          query_sequence='AAAAAA')
        paired_alignment = samtools.PairedAlignment(left, right, tag_length=1)
        self.assertEqual((101,251), paired_alignment.positions)

    def test_filter_value(self):
        left = ConnorAlign(mock_align(), filter_value=None)
        right = ConnorAlign(mock_align(), filter_value=None)
        paired_alignment = samtools.PairedAlignment(left, right, tag_length=1)
        self.assertEqual(None, paired_alignment.filter_value)

        left = ConnorAlign(mock_align(), filter_value='')
        right = ConnorAlign(mock_align(), filter_value='')
        paired_alignment = samtools.PairedAlignment(left, right, tag_length=1)
        self.assertEqual(None, paired_alignment.filter_value)

        left = ConnorAlign(mock_align(), filter_value='foo')
        right = ConnorAlign(mock_align(), filter_value=None)
        paired_alignment = samtools.PairedAlignment(left, right, tag_length=1)
        self.assertEqual(('foo', None), paired_alignment.filter_value)

        left = ConnorAlign(mock_align(), filter_value=None)
        right = ConnorAlign(mock_align(), filter_value='bar')
        paired_alignment = samtools.PairedAlignment(left, right, tag_length=1)
        self.assertEqual((None, 'bar'), paired_alignment.filter_value)

    def test_query_name(self):
        left = mock_align(query_name="alignA", reference_start=100)
        right = mock_align(query_name="alignA", reference_start=200)
        paired_alignment = samtools.PairedAlignment(left, right, tag_length=1)
        self.assertEqual("alignA", paired_alignment.query_name)

    def test_eq(self):
        left = align_seg("alignA", 'chr1', 100, 200, "AAANNNNNNN")
        right = align_seg("alignA", 'chr1', 200, 100, "CCCNNNNNNN")
        other = align_seg("alignA", 'chr2', 100, 200, "AAANNNNNNN")

        base = samtools.PairedAlignment(left, right)
        self.assertEquals(base, samtools.PairedAlignment(left, right))
        self.assertNotEquals(base, samtools.PairedAlignment(other, right))
        self.assertNotEquals(base, samtools.PairedAlignment(left, other))

    def test_hash(self):
        left_A = MockAlignSegment("alignA", 'chr1', 100, 200, "AAATTT" "NNNN")
        right_A = MockAlignSegment("alignA", 'chr1', 200, 100, "NNNN" "CCCGGG")
        left_B = MockAlignSegment("alignA", 'chr1', 100, 200, "AAATTT" "NNNN")
        right_B = MockAlignSegment("alignA", 'chr1', 200, 100, "NNNN" "CCCGGG")

        actual_set = set()
        base = samtools.PairedAlignment(left_A, right_A)
        actual_set.add(base)
        self.assertEquals(1, len(actual_set))

        actual_set.add(base)
        self.assertEquals(1, len(actual_set))

        actual_set.add(samtools.PairedAlignment(left_A, right_A))
        self.assertEquals(1, len(actual_set))

        equivalent_pair = samtools.PairedAlignment(left_B, right_B)
        actual_set.add(equivalent_pair)
        self.assertEquals(1, len(actual_set))

    def test_replace_umt(self):
        # pysam's represtation of the sequence is inconsistent across pysam
        # and python versions; this hack makes the values comparable
        def _barray_to_string(sequence):
            if isinstance(sequence, str):
                return sequence
            else:
                return str(sequence.decode("utf-8"))

        left_A = mock_align(query_sequence='AANN', query_qualities=[1,2,3,4])
        right_A = mock_align(query_sequence='NNCC', query_qualities=[5,6,7,8])
        paired_align = samtools.PairedAlignment(left_A, right_A, tag_length=2)

        paired_align.replace_umt(('GG','TT'))

        self.assertEquals('GGNN',
                          _barray_to_string(paired_align.left.query_sequence))
        self.assertEquals('NNTT',
                          _barray_to_string(paired_align.right.query_sequence))
        self.assertEquals([1,2,3,4],
                          paired_align.left.query_qualities)
        self.assertEquals([5,6,7,8],
                          paired_align.right.query_qualities)

    def test_replace_umt_errorIfInconsistentUmtLength(self):
        left_A = mock_align(query_sequence='AANN', query_qualities=[1,2,3,4])
        right_A = mock_align(query_sequence='NNCC', query_qualities=[5,6,7,8])
        paired_align = samtools.PairedAlignment(left_A, right_A, tag_length=2)

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

    @staticmethod
    def connor_align(query_name, query_sequence, mapping_quality):
        return ConnorAlign(mock_align(query_name=query_name,
                                      query_sequence=query_sequence,
                                      mapping_quality=mapping_quality))

    def test_consensus_qualities_maxMappingQualityScores(self):
        alignAL = self.connor_align('alignA', 'nGT', 30)
        alignAR = self.connor_align('alignA', 'nCT', 25)
        pairA = samtools.PairedAlignment(alignAL, alignAR, tag_length=1)

        alignBL = self.connor_align('alignB', 'nGT', 20)
        alignBR = self.connor_align('alignB', 'nCT', 15)
        pairB = samtools.PairedAlignment(alignBL, alignBR, tag_length=1)

        alignCL = self.connor_align('alignC', 'nGT', 10)
        alignCR = self.connor_align('alignC', 'nCT', 5)
        pairC = samtools.PairedAlignment(alignCL, alignCR, tag_length=1)

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
        pairA = samtools.PairedAlignment(alignAL, alignAR, tag_length=1)

        alignBL = self.connor_align('alignB', "nGT", 30)
        alignBR = self.connor_align('alignB', "nCT", 25)
        pairB = samtools.PairedAlignment(alignBL, alignBR, tag_length=1)

        alignCL = self.connor_align('alignC', "nGT", 10)
        alignCR = self.connor_align('alignC', "nCT", 5)
        pairC = samtools.PairedAlignment(alignCL, alignCR, tag_length=1)

        alignment_pairs = [pairA, pairB, pairC]

        actual_template = connor._TagFamily._select_template_alignment_pair(alignment_pairs)

        self.assertEquals(pairB, actual_template)

    def test_select_template_alignment_pair_breaksTiesByQueryName(self):
        alignAL = self.connor_align('alignA', "nGT", 20)
        alignAR = self.connor_align('alignA', "nCT", 15)
        pairA = samtools.PairedAlignment(alignAL, alignAR, tag_length=1)

        alignBL = self.connor_align('alignB', "nGT", 20)
        alignBR = self.connor_align('alignB', "nCT", 15)
        pairB = samtools.PairedAlignment(alignBL, alignBR, tag_length=1)
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


    def test_init_addsFilterForMinority(self):
        def pair(name, c1, c2):
            left = ConnorAlign(mock_align(query_name=name, cigarstring=c1))
            right = ConnorAlign(mock_align(query_name=name, cigarstring=c2))
            return samtools.PairedAlignment(left, right)
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
            left = ConnorAlign(mock_align(cigarstring=c1))
            right = ConnorAlign(mock_align(cigarstring=c2))
            return samtools.PairedAlignment(left, right)
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


class CoordinateFamilyHolder(BaseConnorTestCase):
    #pylint: disable=no-self-use
    def _pair(self, name, left_pos, right_pos, right_end, reference_name='1'):
        alignL = ConnorAlign(MicroMock(query_name=name,
                                       reference_start=left_pos,
                                       next_reference_start=right_pos,
                                       query_sequence='AAATTTT',
                                       reference_end=0,
                                       reference_name=reference_name))
        alignR = ConnorAlign(MicroMock(query_name=name,
                                       reference_start=right_pos,
                                       next_reference_start=left_pos,
                                       query_sequence='AAATTTT',
                                       reference_end=right_end,
                                       reference_name=reference_name))
        return samtools.PairedAlignment(alignL, alignR)

    def test_build_coordinate_families_noFamilies(self):
        holder = _CoordinateFamilyHolder()
        pairs = []
        actual_families = [f for f in holder.build_coordinate_families(pairs)]
        self.assertEquals([], actual_families)

    def test_build_coordinate_families_gathersSimilarCoordinates(self):
        pair1 = self._pair('1', 100, 200, 205)
        pair2 = self._pair('2', 100, 200, 205)
        pair3 = self._pair('3', 100, 300, 305)
        pair4 = self._pair('4', 100, 300, 305)
        pairs = [pair1, pair2, pair3, pair4]
        holder = _CoordinateFamilyHolder()
        actual_families = set()
        for family in holder.build_coordinate_families(pairs):
            actual_families.add(tuple(family))
        expected_families = set([(pair1, pair2), (pair3, pair4)])
        self.assertEquals(expected_families, actual_families)

    def test_build_coordinate_families_flushesRemaining(self):
        pair1 = self._pair('1', 100, 200, 205)
        pair2 = self._pair('2', 100, 300, 305)
        pair3 = self._pair('3', 200, 400, 405)
        pairs = [pair1, pair3, pair2]
        holder = _CoordinateFamilyHolder()
        actual_families = set()
        for family in holder.build_coordinate_families(pairs):
            actual_families.add(tuple(family))
        expected_families = set([(pair1,), (pair2,), (pair3,)])
        self.assertEquals(expected_families, actual_families)

    def test_build_coordinate_families_streamsPartialResults(self):
        pair1 = self._pair('1', 100, 200, 205)
        pair2 = self._pair('2', 100, 300, 305)

        class AngryError(ValueError):
            pass

        class AngryPair(object):
            @property
            def right(self):
                raise AngryError("I'm angry")
        pair3 = AngryPair()

        pairs = [pair1, pair2, pair3]
        holder = _CoordinateFamilyHolder()
        actual_families = set()
        try:
            for family in holder.build_coordinate_families(pairs):
                actual_families.add(tuple(family))
            self.fail("Expected AngryError")
        except AngryError:
            pass
        expected_families = set([(pair1,),])
        self.assertEquals(expected_families, actual_families)

    def test_pending_pair_count_and_peak(self):
        pair1 = self._pair('A1', 100, 150, 155, reference_name='1')
        pair2 = self._pair('B1', 200, 250, 255, reference_name='1')
        pair3 = self._pair('B2', 200, 250, 355, reference_name='1')

        holder = _CoordinateFamilyHolder()
        self.assertEqual(0, holder.pending_pair_count)
        self.assertEqual(0, holder.pending_pair_peak_count)
        holder._add(pair1)
        self.assertEqual(1, holder.pending_pair_count)
        self.assertEqual(1, holder.pending_pair_peak_count)
        holder._add(pair2)
        self.assertEqual(2, holder.pending_pair_count)
        self.assertEqual(2, holder.pending_pair_peak_count)
        holder._add(pair3)
        self.assertEqual(3, holder.pending_pair_count)
        self.assertEqual(3, holder.pending_pair_peak_count)
        for _ in holder._completed_families('1', 200):
            pass
        self.assertEqual(2, holder.pending_pair_count)
        self.assertEqual(3, holder.pending_pair_peak_count)
        for _ in holder._remaining_families():
            pass
        self.assertEqual(0, holder.pending_pair_count)
        self.assertEqual(3, holder.pending_pair_peak_count)

    def test_build_coordinate_families_6families(self):
        pair1 = self._pair('1', 100, 400, 405)
        pair2 = self._pair('2', 100, 400, 405)
        pair3 = self._pair('3', 200, 400, 405)
        pair4 = self._pair('4', 200, 400, 405)
        pair5 = self._pair('5', 200, 500, 505)
        pair6 = self._pair('6', 300, 600, 605)
        pairs = [pair1, pair2, pair3, pair4, pair5, pair6]

        holder = _CoordinateFamilyHolder()
        actual_coord_families = set()
        for family in holder.build_coordinate_families(pairs):
            actual_coord_families.add(tuple(family))
        expected_coord_families = set([(pair1, pair2),
                                       (pair3, pair4),
                                       (pair5,),
                                       (pair6,)])
        self.assertEqual(expected_coord_families, actual_coord_families)

    def test_build_coordinate_families_flushesInProgressForNewChromosome(self):
        pair1A = self._pair('1A', 100, 200, 205, reference_name='1')
        pair1B = self._pair('1B', 100, 200, 205, reference_name='1')
        pair1C = self._pair('1C', 300, 400, 405, reference_name='1')
        pair2A = self._pair('2A', 100, 200, 205, reference_name='2')
        pair2B = self._pair('2B', 100, 200, 205, reference_name='2')
        pair2C = self._pair('2C', 400, 500, 505, reference_name='2')
        pair3A = self._pair('3A', 600, 700, 705, reference_name='3')
        pair3B = self._pair('3B', 600, 700, 705, reference_name='3')
        pair3C = self._pair('3C', 800, 900, 905, reference_name='3')
        pairs = [pair1A, pair1B, pair1C,
                 pair2A, pair2B, pair2C,
                 pair3A, pair3B, pair3C]
        holder = _CoordinateFamilyHolder()
        actual_families = set()
        for family in holder.build_coordinate_families(pairs):
            actual_families.add(tuple(family))
        expected_families = set([(pair1A, pair1B),
                                 (pair1C,),
                                 (pair2A, pair2B),
                                 (pair2C,),
                                 (pair3A, pair3B),
                                 (pair3C,)])
        self.assertEquals(expected_families, actual_families)


class ConnorTest(BaseConnorTestCase):
    def test_build_coordinate_pairs_singlePair(self):
        align1L = ConnorAlign(mock_align(query_name='1',
                                         reference_start=100,
                                         next_reference_start=200))
        align1R = ConnorAlign(mock_align(query_name='1',
                                         reference_start=200,
                                         next_reference_start=100))
        aligns = [align1L, align1R]

        actual_pairs = [p for p in connor._build_coordinate_pairs(aligns, None)]

        self.assertEqual(1, len(actual_pairs))
        actual_pair = actual_pairs[0]
        self.assertEqual(align1L, actual_pair.left)
        self.assertEqual(align1R, actual_pair.right)

    def test_build_coordinate_pairs_oneCoordinate(self):
        align1L = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=200))
        align2L = ConnorAlign(mock_align(query_name = '2',
                                         reference_start=100,
                                         next_reference_start=200))
        align1R = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=200,
                                         next_reference_start=100))
        align2R = ConnorAlign(mock_align(query_name = '2',
                                         reference_start=200,
                                         next_reference_start=100))
        aligns = [align1L, align2L, align1R, align2R]

        actual_pairs = [f for f in connor._build_coordinate_pairs(aligns, None)]

        self.assertEqual(2, len(actual_pairs))
        pairs = dict([(pair.query_name, pair) for pair in actual_pairs])
        self.assertEqual(align1L, pairs['1'].left)
        self.assertEqual(align1R, pairs['1'].right)
        self.assertEqual(align2L, pairs['2'].left)
        self.assertEqual(align2R, pairs['2'].right)

    def test_build_coordinate_pairs_twoCoordinatesSameRight(self):
        align1L = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=200))
        align2L = ConnorAlign(mock_align(query_name = '2',
                                         reference_start=100,
                                         next_reference_start=200))
        align3L = ConnorAlign(mock_align(query_name = '3',
                                         reference_start=125,
                                         next_reference_start=200))
        align1R = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=200,
                                         next_reference_start=100))
        align2R = ConnorAlign(mock_align(query_name = '2',
                                         reference_start=200,
                                         next_reference_start=100))
        align3R = ConnorAlign(mock_align(query_name = '3',
                                         reference_start=200,
                                         next_reference_start=125))
        aligns = [align1L, align2L, align3L, align1R, align2R, align3R]

        actual_pairs = [f for f in connor._build_coordinate_pairs(aligns, None)]

        self.assertEqual(3, len(actual_pairs))
        actual_pair_names = set([pair.query_name for pair in actual_pairs])
        self.assertEqual(set(['1', '2', '3']), actual_pair_names)

    def test_build_coordinate_pairs_identicalCoordinates(self):
        align1L = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=100))
        align1R = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=100))
        aligns = [align1L, align1R]
        writer = MockAlignWriter()
        actual_pairs = [f for f in connor._build_coordinate_pairs(aligns, writer)]

        self.assertEqual(1, len(actual_pairs))
        actual_pair_names = set([pair.query_name for pair in actual_pairs])
        self.assertEqual(set(['1']), actual_pair_names)


    def test_build_coordinate_pairs_orphanedRightIsSafe(self):
        align1L = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=200))
        align1R = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=200,
                                         next_reference_start=100))
        align2R = ConnorAlign(mock_align(query_name = '2',
                                         reference_start=200,
                                         next_reference_start=100))
        aligns = [align1L, align1R, align2R]
        writer = MockAlignWriter()
        actual_pairs = [f for f in connor._build_coordinate_pairs(aligns, writer)]

        self.assertEqual(1, len(actual_pairs))
        actual_pair_names = set([pair.query_name for pair in actual_pairs])
        self.assertEqual(set(['1']), actual_pair_names)

    def test_build_coordinate_pairs_orphanedRightWrittenToExcluded(self):
        align1L = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=200))
        align2R = ConnorAlign(mock_align(query_name = '2',
                                         reference_start=150,
                                         next_reference_start=100))
        align1R = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=200,
                                         next_reference_start=100))
        aligns = [align1L, align1R, align2R]
        writer = MockAlignWriter()

        for _ in connor._build_coordinate_pairs(aligns, writer):
            pass

        self.assertEqual(1, len(writer._write_calls))
        (actual_family, actual_align) = writer._write_calls[0]
        self.assertEqual(None, actual_family)
        self.assertEqual(align2R, actual_align)
        self.assertEqual('read mate was missing or excluded',
                         actual_align.filter_value)

    def test_build_coordinate_pairs_whenExhasutedRemaindersWrittenToExcluded(self):
        align1L = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=200))
        align1R = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=200,
                                         next_reference_start=100))
        align3L = ConnorAlign(mock_align(query_name = '3',
                                         reference_start=150,
                                         next_reference_start=250))
        aligns = [align1L, align3L, align1R]
        writer = MockAlignWriter()

        for _ in connor._build_coordinate_pairs(aligns, writer):
            pass

        self.assertEqual(1, len(writer._write_calls))
        (actual_family, actual_align) = writer._write_calls[0]
        self.assertEqual(None, actual_family)
        self.assertEqual(align3L, actual_align)
        self.assertEqual('read mate was missing or excluded',
                         actual_align.filter_value)

    def test_build_coordinate_pairs_lookForPassedPops(self):
        align1L = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=200))
        align2L = ConnorAlign(mock_align(query_name = '2',
                                         reference_start=100,
                                         next_reference_start=200))
        align3L = ConnorAlign(mock_align(query_name = '3',
                                         reference_start=300,
                                         next_reference_start=400))
        align1R = ConnorAlign(mock_align(query_name = '1',
                                         reference_start=200,
                                         next_reference_start=100))
        align2R = ConnorAlign(mock_align(query_name = '2',
                                         reference_start=200,
                                         next_reference_start=100))

        class BombError(ValueError):
            pass

        class BombAlign(object):
            @property
            def orientation(self):
                raise BombError("Boom.")
        align4L = BombAlign()

        aligns = [align1L, align2L, align1R, align2R, align3L, align4L]
        actual_pairs = []
        try:
            for pair in connor._build_coordinate_pairs(aligns, None):
                actual_pairs.append(pair)
            self.fail("expected BombError")
        except BombError:
            pass

        self.assertEqual(2, len(actual_pairs))
        actual_pair_names = set([pair.query_name for pair in actual_pairs])
        self.assertEqual(set(['1', '2']), actual_pair_names)


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

    def test_progress_logger_passesThroughItems(self):
        base_gen = [0,1,2,3,4,5,6,7,8,9]
        total_rows = len(base_gen)
        actual_items = []
        for item in connor._progress_logger(base_gen,
                                         total_rows,
                                         self.mock_logger):
            actual_items.append(item)
        self.assertEqual(actual_items, base_gen)

    def test_progress_logger_logsProgress(self):
        base_gen = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
        total_rows = len(base_gen)
        gen = connor._progress_logger(base_gen,
                                     total_rows,
                                     self.mock_logger)
        next(gen)
        self.assertEqual("0% (1/20) alignments processed",
                         self.mock_logger._log_calls['INFO'][0])
        next(gen)
        self.assertEqual("10% (2/20) alignments processed",
                         self.mock_logger._log_calls['INFO'][1])
        next(gen)
        next(gen)
        self.assertEqual("20% (4/20) alignments processed",
                         self.mock_logger._log_calls['INFO'][2])
        next(gen)
        next(gen)
        self.assertEqual("30% (6/20) alignments processed",
                         self.mock_logger._log_calls['INFO'][3])
        for _ in range(7, 21):
            next(gen)

        self.assertRaises(StopIteration,
                          next,
                          gen)

        self.assertEqual("100% (20/20) alignments processed",
                         self.mock_logger._log_calls['INFO'][-1])

    def test_progress_logger_logsMem(self):
        base_gen = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
        total_rows = len(base_gen)
        def supplemental_log(logger):
            logger.debug("foo")
            logger.debug("bar")
        gen = connor._progress_logger(base_gen,
                                     total_rows,
                                     self.mock_logger,
                                     supplemental_log)
        next(gen)
        next(gen)
        self.assertRegexpMatches(self.mock_logger._log_calls['DEBUG'][0],
                                 "foo")
        self.assertRegexpMatches(self.mock_logger._log_calls['DEBUG'][1],
                                 "bar")
        next(gen)
        self.assertRegexpMatches(self.mock_logger._log_calls['DEBUG'][2],
                                 "foo")
        self.assertRegexpMatches(self.mock_logger._log_calls['DEBUG'][3],
                                 "bar")
        next(gen)
        self.assertRegexpMatches(self.mock_logger._log_calls['DEBUG'][4],
                                 "foo")
        self.assertRegexpMatches(self.mock_logger._log_calls['DEBUG'][5],
                                 "bar")



    def test_rank_tags_breaksTiesByTag(self):
        pair0 = align_pair("align0", 'chr1', 100, 200, "TTTNNN", "NNNGGG")
        pair1 = align_pair("align1", 'chr1', 100, 200, "AAANNN", "NNNCCC")
        pair2 = align_pair("align2", 'chr1', 100, 200, "AAANNN", "NNNGGG")
        input_aligns = [pair0, pair1, pair2]

        actual_tags = connor._rank_tags(input_aligns)

        expected_tags = [('AAA', 'CCC'), ('AAA', 'GGG'), ('TTT', 'GGG')]
        self.assertEquals(expected_tags, actual_tags)


class ConnorIntegrationTestCase(BaseConnorTestCase):
    def setUp(self):
        self.check_sysout_safe()
        BaseConnorTestCase.setUp(self)

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
