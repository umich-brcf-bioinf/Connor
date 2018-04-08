#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments
#pylint: disable=deprecated-method

from __future__ import print_function, absolute_import, division

from test.utils_test import BaseConnorTestCase
from test.utils_test import MicroMock

from connor.family import CoordinateFamilyHolder
from connor.consam.alignments import ConnorAlign
from connor.consam.alignments import PairedAlignment
from connor.family import TagFamily


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

# #TODO: cgates: replace this with self.mock_align
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

class TagFamiliyTest(BaseConnorTestCase):
    def align_pair(self, q, rn, rs, nrs, s1, s2, tag_length=3):
        alignL = self.align_seg(q, rn, rs, nrs, s1)
        alignR = self.align_seg(q, rn, rs, nrs, s2)
        return PairedAlignment(alignL, alignR, tag_length)

    @staticmethod
    def align_seg( #pylint: disable=dangerous-default-value
                  query_name,
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


    def test_init(self):
        pair1 = self.align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        pair2 = self.align_pair('alignB', 'chr1', 100, 200, 'GGGNNN', 'TTTNNN')
        pair3 = self.align_pair('alignC', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        align_pairs = [pair1, pair2, pair3]
        inexact_match_count = 42

        input_umt = "AAANNNCCCNNN"
        actual_tag_family = TagFamily(input_umt,
                                             align_pairs,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(input_umt, actual_tag_family.umt())
        self.assertEquals(align_pairs, actual_tag_family.align_pairs)
        self.assertEquals(42, actual_tag_family.inexact_match_count)

    def test_init_setsFilterValue(self):
        pair1 = self.align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        pair2 = self.align_pair('alignB', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        pair3 = self.align_pair('alignC', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        align_pairs = [pair1, pair2, pair3]
        inexact_match_count = 42
        input_umt = "AAANNNCCCNNN"

        def family_filter(family):
            return ":".join([a.query_name for a in family.align_pairs]) +\
                 " is too odd"
        actual_tag_family = TagFamily(input_umt,
                                              align_pairs,
                                              inexact_match_count,
                                              consensus_threshold=0.6,
                                              family_filter=family_filter)

        self.assertEquals('alignA:alignB:alignC is too odd',
                          actual_tag_family.filter_value)

    def test_is_consensus_template(self):
        def pair(name):
            left = ConnorAlign(self.mock_align(query_name=name))
            right = ConnorAlign(self.mock_align(query_name=name))
            return PairedAlignment(left, right)
        pairA = pair('alignA')
        alignA = ConnorAlign(self.mock_align(query_name='alignA'))
        alignB = ConnorAlign(self.mock_align(query_name='alignB'))
        alignC = ConnorAlign(self.mock_align(query_name='alignC'))
        tag_family = TagFamily(umt=("AAA", "CCC"),
                                       alignments=[pairA],
                                       inexact_match_count=0,
                                       consensus_threshold=0.6)

        self.assertEqual(True, tag_family.is_consensus_template(alignA))
        self.assertEqual(False, tag_family.is_consensus_template(alignB))
        self.assertEqual(False, tag_family.is_consensus_template(alignC))

    def test_umt(self):
        left = ConnorAlign(self.mock_align())
        right = ConnorAlign(self.mock_align())
        pair1 = PairedAlignment(left, right)

        input_umt = ('AAANNN', 'CCCNNN')
        actual_tag_family = TagFamily(input_umt,
                                             [pair1],
                                             inexact_match_count=42,
                                             consensus_threshold=0.6)

        self.assertEquals(input_umt, actual_tag_family.umt())
        self.assertEquals('AAANNN~CCCNNN',
                          actual_tag_family.umt(format_string="{left}~{right}"))

    def test_consensus_rewrites_umt(self):
        pair1 = self.align_pair('alignA', 'chr1', 100, 200, 'GGGNNN', 'NNNTTT')
        input_umts = ("AAA", "CCC")
        inexact_match_count = 0

        actual_tag_family = TagFamily(input_umts,
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
        pair1 = self.align_pair('alignA', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        pair2 = self.align_pair('alignB', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        pair3 = self.align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)
        actual_consensus_seq = actual_tag_family.consensus

        self.assertEquals("nnnGGG",
                          actual_consensus_seq.left.query_sequence)
        self.assertEquals("TTTnnn",
                          actual_consensus_seq.right.query_sequence)

    def test_consensus_sequence_majority_wins(self):
        pair1 = self.align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'TCTnnn')
        pair2 = self.align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'TCTnnn')
        pair3 = self.align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")

        actual_tag_family = TagFamily(input_umts,
                                              alignments,
                                              inexact_match_count=1,
                                              consensus_threshold=0.5)
        consensus_pair = actual_tag_family.consensus

        self.assertEquals("nnnGTG",
                          consensus_pair.left.query_sequence)
        self.assertEquals("TCTnnn",
                          consensus_pair.right.query_sequence)

    def test_consensus_sequence_below_threshold_Ns(self):
        pair1 = self.align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'TCTnnn')
        pair2 = self.align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'TCTnnn')
        pair3 = self.align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'TTTnnn')
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        threshold = 0.8
        inexact_match_count = 0

        actual_tag_family = TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             threshold)
        consensus_pair = actual_tag_family.consensus

        self.assertEquals("nnnGNG",
                          consensus_pair.left.query_sequence)
        self.assertEquals("TNTnnn",
                          consensus_pair.right.query_sequence)

    def connor_align(self, query_name, query_sequence, mapping_quality):
        return ConnorAlign(self.mock_align(query_name=query_name,
                                      query_sequence=query_sequence,
                                      mapping_quality=mapping_quality))

    def test_consensus_qualities_maxMappingQualityScores(self):
        alignAL = self.connor_align('alignA', 'nGT', 30)
        alignAR = self.connor_align('alignA', 'nCT', 25)
        pairA = PairedAlignment(alignAL, alignAR, tag_length=1)

        alignBL = self.connor_align('alignB', 'nGT', 20)
        alignBR = self.connor_align('alignB', 'nCT', 15)
        pairB = PairedAlignment(alignBL, alignBR, tag_length=1)

        alignCL = self.connor_align('alignC', 'nGT', 10)
        alignCR = self.connor_align('alignC', 'nCT', 5)
        pairC = PairedAlignment(alignCL, alignCR, tag_length=1)

        alignments = [pairA, pairB, pairC]
        input_umts = ("n", "n")
        inexact_match_count = 0

        actual_tag_family = TagFamily(input_umts,
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
        pairA = PairedAlignment(alignAL, alignAR, tag_length=1)

        alignBL = self.connor_align('alignB', "nGT", 30)
        alignBR = self.connor_align('alignB', "nCT", 25)
        pairB = PairedAlignment(alignBL, alignBR, tag_length=1)

        alignCL = self.connor_align('alignC', "nGT", 10)
        alignCR = self.connor_align('alignC', "nCT", 5)
        pairC = PairedAlignment(alignCL, alignCR, tag_length=1)

        alignment_pairs = [pairA, pairB, pairC]

        actual_template = TagFamily._select_template_alignment_pair(alignment_pairs)

        self.assertEquals(pairB, actual_template)

    def test_select_template_alignment_pair_breaksTiesByQueryName(self):
        alignAL = self.connor_align('alignA', "nGT", 20)
        alignAR = self.connor_align('alignA', "nCT", 15)
        pairA = PairedAlignment(alignAL, alignAR, tag_length=1)

        alignBL = self.connor_align('alignB', "nGT", 20)
        alignBR = self.connor_align('alignB', "nCT", 15)
        pairB = PairedAlignment(alignBL, alignBR, tag_length=1)
        alignment_pairs = [pairA, pairB]

        actual_template = TagFamily._select_template_alignment_pair(alignment_pairs)

        self.assertEquals(pairA, actual_template)


    def test_consensus_uniform_cigars_admitted(self):
        pair1 = self.align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left.cigarstring = "3S3M"
        pair1.right.cigarstring = "3S3M"
        pair2 = self.align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left.cigarstring = "3S3M"
        pair2.right.cigarstring = "3S3M"
        pair3 = self.align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left.cigarstring = "3S3M"
        pair3.right.cigarstring = "3S3M"
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        included_names = [a.left.query_name for a in actual_tag_family.align_pairs]
        self.assertEquals(['alignA','alignB','alignC'], included_names)

    def test_consensus_minority_cigars_excluded(self):
        pair1 = self.align_pair('alignA', 1, 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left.cigarstring = "3S3M"
        pair1.right.cigarstring = "3S3M"
        pair2 = self.align_pair('alignB', 1, 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left.cigarstring = "3S3M"
        pair2.right.cigarstring = "3S3M"
        pair3 = self.align_pair('alignC', 1, 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left.cigarstring = "3S3I"
        pair3.right.cigarstring = "3S3M"
        align_pairs = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = TagFamily(input_umts,
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
        pair1 = self.align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left.cigarstring = "3S3M"
        pair1.right.cigarstring = "3S3M"
        pair2 = self.align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left.cigarstring = "3S3M"
        pair2.right.cigarstring = "3S3M"
        pair3 = self.align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left.cigarstring = "3S3I"
        pair3.right.cigarstring = "3S3I"
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(2, actual_tag_family.distinct_cigar_count)

    def test_init_distinctCigarCountOneForConsisentCigar(self):
        pair1 = self.align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left.cigarstring = "3S3M"
        pair1.right.cigarstring = "3S3M"
        pair2 = self.align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left.cigarstring = "3S3M"
        pair2.right.cigarstring = "3S3M"
        pair3 = self.align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left.cigarstring = "3S3M"
        pair3.right.cigarstring = "3S3M"
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(1, actual_tag_family.distinct_cigar_count)

    def test_distinct_cigar_count_single_ended(self):
        pair1 = self.align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair1.left.cigarstring = "3S3M"
        pair1.right.cigarstring = "3S3M"
        pair2 = self.align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2.left.cigarstring = "3S3M"
        pair2.right.cigarstring = "3S3M"
        pair3 = self.align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3.left.cigarstring = "3S3I"
        pair3.right.cigarstring = "3S3M"
        alignments = [pair1, pair2, pair3]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(2, actual_tag_family.distinct_cigar_count)


    def test_init_addsFilterForMinority(self):
        def pair(name, c1, c2):
            left = ConnorAlign(self.mock_align(query_name=name, cigarstring=c1))
            right = ConnorAlign(self.mock_align(query_name=name, cigarstring=c2))
            return PairedAlignment(left, right)
        pairA = pair('alignA', '3S3M', '3S3M')
        pairB = pair('alignB', '3S3M', '3S3M')
        pairC = pair('alignC', '3S1I3M', '3S3M')
        align_pairs = [pairA, pairB, pairC]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = TagFamily(input_umts,
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
            left = ConnorAlign(self.mock_align(cigarstring=c1))
            right = ConnorAlign(self.mock_align(cigarstring=c2))
            return PairedAlignment(left, right)
        pairA = pair("3S3M", "3S3M")
        pairB = pair("3S3M", "3S3M")
        pairC = pair("3S1I3M", "3S3M")
        alignments = [pairA, pairB, pairC]
        input_umts = ("nnn", "nnn")
        inexact_match_count = 0

        actual_tag_family = TagFamily(input_umts,
                                             alignments,
                                             inexact_match_count,
                                             consensus_threshold=0.6)

        self.assertEquals(2, actual_tag_family.included_pair_count)


class CoordinateFamilyHolderTest(BaseConnorTestCase):
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
        return PairedAlignment(alignL, alignR)

    def test_build_coordinate_families_noFamilies(self):
        holder = CoordinateFamilyHolder()
        pairs = []
        actual_families = [f for f in holder.build_coordinate_families(pairs)]
        self.assertEquals([], actual_families)

    def test_build_coordinate_families_gathersSimilarCoordinates(self):
        pair1 = self._pair('1', 100, 200, 205)
        pair2 = self._pair('2', 100, 200, 205)
        pair3 = self._pair('3', 100, 300, 305)
        pair4 = self._pair('4', 100, 300, 305)
        pairs = [pair1, pair2, pair3, pair4]
        holder = CoordinateFamilyHolder()
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
        holder = CoordinateFamilyHolder()
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
        holder = CoordinateFamilyHolder()
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

        holder = CoordinateFamilyHolder()
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

        holder = CoordinateFamilyHolder()
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
        holder = CoordinateFamilyHolder()
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
