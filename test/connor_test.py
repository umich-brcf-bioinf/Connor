#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments
from __future__ import print_function, absolute_import, division
from collections import namedtuple
import os
import unittest
import pysam
from testfixtures.tempdirectory import TempDirectory
from connor import connor
from connor import samtools

class MockLogger(object):
    def __init__(self):
        self._log_calls = []

    def log(self, msg_format, *args):
        self._log_calls.append((msg_format, args))

class MicroMock(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):
        return 1

class MockAlignSegment(object):
    def __init__(self,
                 query_name,
                 reference_name,
                 reference_start,
                 next_reference_start,
                 query_sequence='AAACCC',
                 query_qualities=[25,25,25,25,25,25],
                 cigarstring='6M',
                 reference_end=None):
        if not reference_end:
            reference_end = reference_start + len(query_sequence)
        self.query_name = query_name
        self.reference_name = reference_name
        self.reference_start = reference_start
        self.next_reference_start = next_reference_start
        self.query_sequence = query_sequence
        self.query_qualities = query_qualities
        self.cigarstring = cigarstring
        self.reference_end = reference_end

    def __hash__(self):
        return hash(self.query_name)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def set_tag(self, name, value, tag_type):
        pass

#TODO: cgates: inline this method
def align_seg(query_name,
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
    
#     if not reference_end:
#         reference_end = reference_start + len(query_sequence)
#     return MicroMock(query_name=query_name,
#                      reference_name=reference_name,
#                      reference_start=reference_start,
#                      next_reference_start=next_reference_start,
#                      query_sequence=query_sequence,
#                      query_qualities=query_qualities,
#                      cigarstring=cigarstring,
#                      reference_end=reference_end)

def align_pair(q, rn, rs, nrs, s1, s2, tag_length=3):
    alignL = align_seg(q, rn, rs, nrs, s1)
    alignR = align_seg(q, rn, rs, nrs, s2)
    return connor.PairedAlignment(alignL, alignR, tag_length)

def _create_file(path, filename, contents):
    filename = os.path.join(path, filename)
    with open(filename, 'wt') as new_file:
        new_file.write(contents)
        new_file.flush()
    return filename

def _create_bam(path, filename, sam_contents, index=True):
    sam_filename = _create_file(path, filename, sam_contents)
    bam_filename = sam_filename.replace(".sam", ".bam")
    _pysam_bam_from_sam(sam_filename, bam_filename, index)
    return bam_filename

def _pysam_bam_from_sam(sam_filename, bam_filename, index=True):
    infile = pysam.AlignmentFile(sam_filename, "r")
    outfile = pysam.AlignmentFile(bam_filename, "wb", template=infile)
    for s in infile:
        outfile.write(s)
    infile.close()
    outfile.close()
    if index:
        samtools.index(bam_filename)

def _pysam_alignments_from_bam(bam_filename):
    infile = pysam.AlignmentFile(bam_filename, "rb")
    aligned_segments = [s for s in infile]
    infile.close()
    return aligned_segments

class PairedAlignmentTest(unittest.TestCase):
    def test_init(self):
        left_align = MockAlignSegment("alignA", 'chr1', 10, 20, "AAATTT" "GGGG", reference_end=14)
        right_align = MockAlignSegment("alignA", 'chr1', 20, 10, "TTTT" "CCCGGG" ,reference_end=24)
        tag_length = 6
        actual_paired_alignment = connor.PairedAlignment(left_align,
                                                         right_align,
                                                         tag_length)

        self.assertIs(left_align, actual_paired_alignment.left_alignment)
        self.assertIs(right_align, actual_paired_alignment.right_alignment)
        self.assertEquals(("AAATTT","CCCGGG"), actual_paired_alignment.umi)

    def test_eq(self):
        left = align_seg("alignA", 'chr1', 100, 200, "AAANNNNNNN")
        right = align_seg("alignA", 'chr1', 200, 100, "CCCNNNNNNN")
        other = align_seg("alignB", 'chr1', 100, 200, "AAANNNNNNN")

        base = connor.PairedAlignment(left, right)
        self.assertEquals(base, connor.PairedAlignment(left, right))
        self.assertNotEquals(base, connor.PairedAlignment(other, right))
        self.assertNotEquals(base, connor.PairedAlignment(left, other))

    def test_hash(self):
        left_A = MockAlignSegment("alignA", 'chr1', 100, 200, "AAATTT" "NNNN")
        right_A = MockAlignSegment("alignA", 'chr1', 200, 100, "NNNN" "CCCGGG")
        left_B = MockAlignSegment("alignA", 'chr1', 100, 200, "AAATTT" "NNNN")
        right_B = MockAlignSegment("alignA", 'chr1', 200, 100, "NNNN" "CCCGGG")

        actual_set = set()
        base = connor.PairedAlignment(left_A, right_A)
        actual_set.add(base)
        self.assertEquals(1, len(actual_set))

        actual_set.add(base)
        self.assertEquals(1, len(actual_set))

        actual_set.add(connor.PairedAlignment(left_A, right_A))
        self.assertEquals(1, len(actual_set))

        equivalent_pair = connor.PairedAlignment(left_B, right_B)
        actual_set.add(equivalent_pair)
        self.assertEquals(1, len(actual_set))

    def test_replace_umi(self):
        left_A = align_seg("alignA", 'chr1', 100, 200, "AAANNNNNNN")
        right_A = align_seg("alignA", 'chr1', 200, 100, "CCCNNNNNNN")
        paired_align = connor.PairedAlignment(left_A, right_A, tag_length=4)

        paired_align.replace_umi(('GGGG','TTTT'))

        self.assertEquals('GGGGNNNNNN',
                          paired_align.left_alignment.query_sequence)
        self.assertEquals('TTTTNNNNNN',
                          paired_align.right_alignment.query_sequence)


class TagFamiliyTest(unittest.TestCase):
    
    def test_build_tag_family(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'GGGNNN', 'TTTNNN')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        alignments = [pair1, pair2, pair3]

        input_umi = "AAANNNCCCNNN"
        actual_tag_family = connor.TagFamily(input_umi, alignments)

        self.assertEquals(input_umi, actual_tag_family.umi)
        self.assertEquals(alignments, actual_tag_family.alignments)


    def test_consensus_enforce_tag(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'GGGGGG', 'CCCNNN')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'GGGNNN', 'TTTNNN')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'CCCNNN', 'CCCTTT')
        alignments = [pair1, pair2, pair3]
        input_umis = ("AAA", "CCC")

        actual_tag_family = connor.TagFamily(input_umis, alignments)
        actual_consensus_seq = actual_tag_family.consensus

        self.assertEquals(input_umis, actual_consensus_seq.umi)

    def test_consensus_sequence_trivial_noop(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")

        actual_tag_family = connor.TagFamily(input_umis, alignments)
        actual_consensus_seq = actual_tag_family.consensus

        self.assertEquals("nnnGGG",
                          actual_consensus_seq.left_alignment.query_sequence)
        self.assertEquals("nnnTTT",
                          actual_consensus_seq.right_alignment.query_sequence)

    def test_consensus_sequence_majority_wins(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        threshold = 1 / 2

        actual_tag_family = connor.TagFamily(input_umis, alignments, threshold)
        consensus_pair = actual_tag_family.consensus

        self.assertEquals("nnnGTG",
                          consensus_pair.left_alignment.query_sequence)
        self.assertEquals("nnnTCT",
                          consensus_pair.right_alignment.query_sequence)

    def test_consensus_sequence_below_threshold_Ns(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'nnnGTG', 'nnnTCT')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'nnnGGG', 'nnnTTT')
        alignments = [pair1, pair2, pair3]
        input_umis = ("nnn", "nnn")
        threshold = 0.8

        actual_tag_family = connor.TagFamily(input_umis, alignments, threshold)
        consensus_pair = actual_tag_family.consensus

        self.assertEquals("nnnGNG",
                          consensus_pair.left_alignment.query_sequence)
        self.assertEquals("nnnTNT",
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

        actual_tag_family = connor.TagFamily(input_umis, alignments)
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

        actual_tag_family = connor.TagFamily(input_umis, alignments)

        actual_names = [a.left_alignment.query_name for a in actual_tag_family.alignments]
        self.assertEquals(['alignA','alignB','alignC'], actual_names)

    def test_consensus_minority_cigars_rejected(self):
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

        actual_tag_family = connor.TagFamily(input_umis, alignments)

        actual_names = [a.left_alignment.query_name for a in actual_tag_family.alignments]
        self.assertEquals(['alignA','alignB'], actual_names)


class ConnorTest(unittest.TestCase):

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

        pair_A = connor.PairedAlignment(align_A0, align_A1)
        pair_B = connor.PairedAlignment(align_B0, align_B1)
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

        actual_families = [family for family in connor._build_coordinate_families(alignments, coord_read_name_manifest)]

        pair_A = connor.PairedAlignment(align_A0, align_A1)
        pair_B = connor.PairedAlignment(align_B0, align_B1)
        pair_C = connor.PairedAlignment(align_C0, align_C1)
        pair_D = connor.PairedAlignment(align_D0, align_D1)

        expected_families = [set([pair_A, pair_B]),
                             set([pair_D]),
                             set([pair_C])]
        self.assertEquals(expected_families, actual_families)

#     def test_build_consensus_pair(self):
#         align_A0 = align_seg("alignA", 'chr1', 10, 100)
#         align_A1 = align_seg("alignA", 'chr1', 100, 10)
#         align_B0 = align_seg("alignB", 'chr1', 10, 100)
#         align_B1 = align_seg("alignB", 'chr1', 100, 10)
#         alignments = set([connor.PairedAlignment(align_A0, align_A1),
#                           connor.PairedAlignment(align_B0, align_B1)])
# 
#         actual_pair = connor._build_consensus_pair(alignments)
# 
#         expected_pair = connor.PairedAlignment(align_B0, align_B1)
#         self.assertEquals(expected_pair, actual_pair)

    def test_build_tag_families(self):
        pair1 = align_pair('alignA', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        pair2 = align_pair('alignB', 'chr1', 100, 200, 'GGGNNN', 'TTTNNN')
        pair3 = align_pair('alignC', 'chr1', 100, 200, 'AAANNN', 'CCCNNN')
        paired_aligns = [pair1, pair2, pair3]
        tag1 = ('AAA', 'CCC')
        tag2 = ('GGG', 'TTT')
        ranked_tags = [tag1, tag2]

        actual_tag_fam_list = connor._build_tag_families(paired_aligns,
                                                         ranked_tags)
        actual_tag_fam = {fam.umi:fam for fam in actual_tag_fam_list}

        self.assertEquals(2, len(actual_tag_fam))
        self.assertEquals(set([pair1, pair3]),
                          set(actual_tag_fam[tag1].alignments))
        self.assertEquals(set([pair2]),
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
                                                         ranked_tags)
        actual_tag_fam = {fam.umi:fam for fam in actual_tag_fam_list}

        self.assertEquals(set([pair1, pair2, pair3, pair4]),
                          set(actual_tag_fam[('AAA','GGG')].alignments))
        self.assertEquals(set([pair5]),
                          set(actual_tag_fam[('AAA','CCC')].alignments))

    def test_sort_and_index_bam(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameB1|147|chr10|400|0|5M|=|200|100|CCCCC|>>>>>
readNameA1|147|chr10|300|0|5M|=|100|100|AAAAA|>>>>>
readNameA1|99|chr10|100|0|5M|=|300|200|AAAAA|>>>>>
readNameB1|99|chr10|200|0|5M|=|400|200|CCCCC|>>>>>
readNameA2|147|chr10|300|0|5M|=|100|100|AAAAA|>>>>>
readNameA2|99|chr10|100|0|5M|=|300|200|AAAAA|>>>>>
'''.replace("|", "\t")

        with TempDirectory() as tmp_dir:
            bam = _create_bam(tmp_dir.path, 
                              "input.sam",
                              sam_contents,
                              index=False)
            connor._sort_and_index_bam(bam)
            alignments = pysam.AlignmentFile(bam, "rb").fetch()
            aligns = [(a.query_name, a.reference_start + 1) for a in alignments]
            self.assertEquals(6, len(aligns))
            self.assertEquals([("readNameA1", 100),
                               ("readNameA2", 100),
                               ("readNameB1", 200),
                               ("readNameA1", 300),
                               ("readNameA2", 300),
                               ("readNameB1", 400)],
                              aligns)

            original_dir = os.getcwd()
            try:
                os.chdir(tmp_dir.path)
                os.mkdir("tmp")
                bam = _create_bam(os.path.join(tmp_dir.path, "tmp"),
                                  "input.sam",
                                  sam_contents,
                                  index=False)
                bam_filename = os.path.basename(bam)

                connor._sort_and_index_bam(os.path.join("tmp", bam_filename))

                alignments = pysam.AlignmentFile(bam, "rb").fetch()
                aligns = [(a.query_name, a.reference_start + 1) for a in alignments]
                self.assertEquals(6, len(aligns))
                self.assertEquals([("readNameA1", 100),
                                   ("readNameA2", 100),
                                   ("readNameB1", 200),
                                   ("readNameA1", 300),
                                   ("readNameA2", 300),
                                   ("readNameB1", 400)],
                                  aligns)
            finally:
                os.chdir(original_dir)



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
        self.assertRaises(connor._ConnorUsageError,
                          connor._parse_command_line_args,
                          ["input"])
        self.assertRaises(connor._ConnorUsageError,
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


class TestLightweightAlignment(unittest.TestCase):
    def test_lightweight_alignment_forwardRead(self):
        alignedSegment = align_seg("align1", 'chr1', 10, 100)

        actual_lwa = connor.LightweightAlignment(alignedSegment)

        self.assertEquals("align1", actual_lwa.name)
        self.assertEquals(('chr1', 10, 100), actual_lwa.key)

    def test_lightweight_alignment_reverseRead(self):
        alignedSegment = align_seg("align1", 'chr1', 100, 10)

        actual_lwa = connor.LightweightAlignment(alignedSegment)

        self.assertEquals("align1", actual_lwa.name)
        self.assertEquals(('chr1', 10, 100), actual_lwa.key)

    def test_lightweight_alignment_weirdRead(self):
        alignedSegment = align_seg("align1", 'chr1', 100, 100)

        actual_lwa = connor.LightweightAlignment(alignedSegment)

        self.assertEquals("align1", actual_lwa.name)
        self.assertEquals(('chr1', 100, 100), actual_lwa.key)


class ConnorIntegrationTestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.original_logger = connor._log
        self.mock_logger = MockLogger()
        connor._log = self.mock_logger.log
        self.min_orig_reads =  connor.MIN_ORIG_READS
        connor.MIN_ORIG_READS = 0

    def tearDown(self):
        connor._log = self.original_logger
        connor.MIN_ORIG_READS = self.min_orig_reads
        unittest.TestCase.tearDown(self)

    def test(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|0|5M|=|300|200|AAAAA|>>>>>
readNameA2|99|chr10|100|0|5M|=|300|200|AAAAA|>>>>>
readNameB1|99|chr10|200|0|5M|=|400|200|CCCCC|>>>>>
readNameA1|147|chr10|300|0|5M|=|100|100|AAAAA|>>>>>
readNameA2|147|chr10|300|0|5M|=|100|100|AAAAA|>>>>>
readNameB1|147|chr10|400|0|5M|=|200|100|CCCCC|>>>>>
'''.replace("|", "\t")

        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, "input.sam", sam_contents)
            output_bam = os.path.join(tmp_dir.path, "output.bam")
            connor.main(["connor", input_bam, output_bam])

            alignments = pysam.AlignmentFile(output_bam, "rb").fetch()

            aligns = [(a.query_name, a.reference_start + 1) for a in alignments]
            self.assertEquals(4, len(aligns))
            self.assertEquals([("readNameA1", 100),
                               ("readNameB1", 200),
                               ("readNameA1", 300),
                               ("readNameB1", 400)],
                              aligns)

    def test_logging(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|0|5M|=|300|200|AAAAA|>>>>>
readNameA2|99|chr10|100|0|5M|=|300|200|AAAAA|>>>>>
readNameB1|99|chr10|200|0|5M|=|400|200|CCCCC|>>>>>
readNameA1|147|chr10|300|0|5M|=|100|100|AAAAA|>>>>>
readNameA2|147|chr10|300|0|5M|=|100|100|AAAAA|>>>>>
readNameB1|147|chr10|400|0|5M|=|200|100|CCCCC|>>>>>
'''.replace("|", "\t")

        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, 'input.sam', sam_contents)
            output_bam = os.path.join(tmp_dir.path, 'output.bam')
            connor.main(["connor", input_bam, output_bam])
            log_calls = self.mock_logger._log_calls
            self.assertEquals("connor begins", log_calls[0][0])
            self.assertEquals((input_bam,), log_calls[1][1])
            self.assertRegexpMatches(log_calls[2][0], 'original read count:')
            self.assertEquals((6,), log_calls[2][1])
            self.assertRegexpMatches(log_calls[3][0], 'consensus read count:')
            self.assertEquals((4,), log_calls[3][1])
            self.assertRegexpMatches(log_calls[4][0],
                                     'consensus/original:')
            self.assertEquals((4 / 6,), log_calls[4][1])
            self.assertEquals("sorting and indexing bam", log_calls[5][0])
            self.assertEquals((output_bam,), log_calls[6][1])
            self.assertEquals("connor complete", log_calls[7][0])

    def test_distinctPairStartsAreNotCombined(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|0|5M|=|300|200|AAAAA|>>>>>
readNameB1|99|chr10|100|0|5M|=|500|200|AAAAA|>>>>>
readNameA1|147|chr10|300|0|5M|=|100|200|AAAAA|>>>>>
readNameB1|147|chr10|500|0|5M|=|100|200|AAAAA|>>>>>
'''.replace("|", "\t")

        with TempDirectory() as tmp_dir:
            input_bam = _create_bam(tmp_dir.path, "input.sam", sam_contents)
            output_bam = os.path.join(tmp_dir.path, "output.bam")
            connor.main(["connor", input_bam, output_bam])
            
            alignments = _pysam_alignments_from_bam(output_bam)
            
            aligns = [(a.query_name, a.reference_start + 1) for a in alignments]
            self.assertEquals(4, len(aligns))
            self.assertEquals([("readNameA1", 100),
                               ("readNameB1", 100),
                               ("readNameA1", 300),
                               ("readNameB1", 500)],
                              aligns)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()