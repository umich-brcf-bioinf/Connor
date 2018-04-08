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

from test.consam_test.writers_test import MockAlignWriter
from test.utils_test import BaseConnorTestCase
from test.utils_test import MicroMock
from testfixtures.tempdirectory import TempDirectory

import connor.connor as connor
import connor.consam.pysamwrapper as pysamwrapper
import connor.consam.writers as writers
from connor.consam.alignments import ConnorAlign
from connor.consam.alignments import PairedAlignment


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
    return PairedAlignment(alignL, alignR, tag_length)


class ConnorTest(BaseConnorTestCase):
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
        actual_tag_fam = {fam.umt():fam for fam in actual_tag_fam_list}

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
        actual_tag_fam = {fam.umt():fam for fam in actual_tag_fam_list}

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
            input_bam = self.create_bam(tmp_dir.path,
                                        "input.sam",
                                        sam_contents)
            output_bam = os.path.join(tmp_dir.path, "output.bam")
            args = Namespace(simplify_pg_header=True,
                             original_command_line='foo')
            consensus_writer = writers.build_writer(input_bam,
                                                     output_bam,
                                                     tags=[],
                                                     args=args)
            annotated_writer = MockAlignWriter()
            args = Namespace(input_bam=input_bam,
                             consensus_freq_threshold=0.6,
                             min_family_size_threshold=0,
                             umt_distance_threshold=1)
            connor._dedup_alignments(args,
                                     consensus_writer,
                                     annotated_writer,
                                     self.mock_logger)
            consensus_writer.close()
            alignments = pysamwrapper.alignment_file(output_bam, "rb").fetch()

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
            input_bam = self.create_bam(tmp_dir.path,
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
                                     MockAlignWriter(),
                                     MockAlignWriter(),
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
            input_bam = self.create_bam(tmp_dir.path,
                                        "input.sam",
                                        sam_contents)
            output_bam = os.path.join(tmp_dir.path, "output.bam")
            args = Namespace(input_bam=input_bam,
                             consensus_freq_threshold=0.6,
                             min_family_size_threshold=0,
                             umt_distance_threshold=1,
                             simplify_pg_header=False,
                             original_command_line='foo')
            consensus_writer = writers.build_writer(input_bam,
                                                     output_bam,
                                                     [],
                                                     args)
            annotated_writer = writers.AlignWriter.NULL

            connor._dedup_alignments(args,
                                     consensus_writer,
                                     annotated_writer,
                                     self.mock_logger)
            consensus_writer.close()
            alignments = self.pysam_alignments_from_bam(output_bam)

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
            input_bam = self.create_bam(tmp_dir.path,
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
