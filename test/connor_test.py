#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
from __future__ import print_function, absolute_import
from collections import namedtuple
import os
import pysam
from testfixtures.tempdirectory import TempDirectory
import unittest
from connor import connor

MockPysamAlignedSegment = namedtuple('MockPysamAlignedSegment',
                                     ('query_name,'
                                      'reference_name,'
                                      'reference_start,'
                                      'next_reference_start,'
                                      'sequence'))

def align_seg(a, b, c, d, e='AAACCC'):
    return MockPysamAlignedSegment(query_name=a,
                                   reference_name=b,
                                   reference_start=c,
                                   next_reference_start=d,
                                   sequence=e)

def _create_file(path, filename, contents):
    filename = os.path.join(path, filename)
    with open(filename, 'wt') as new_file:
        new_file.write(contents)
        new_file.flush()
    return filename

def _create_bam(path, filename, sam_contents):
    sam_filename = _create_file(path, filename, sam_contents)
    bam_filename = sam_filename.replace(".sam", ".bam")
    _pysam_bam_from_sam(sam_filename, bam_filename)
    return bam_filename

def _pysam_bam_from_sam(sam_filename, bam_filename):
    infile = pysam.AlignmentFile(sam_filename, "r")
    outfile = pysam.AlignmentFile(bam_filename, "wb", template=infile)
    for s in infile:
        outfile.write(s)
    infile.close()
    outfile.close()
    pysam.index(bam_filename, catch_stdout=False)

def _pysam_alignments_from_bam(bam_filename):
    infile = pysam.AlignmentFile(bam_filename, "rb")
    aligned_segments = [s for s in infile]
    infile.close()
    return aligned_segments

class PairedAlignmentTest(unittest.TestCase):
    def test_init(self):
        alignment1l = align_seg("alignA", 'chr1', 100, 200, "AAANNNNNNN")
        alignment1r = align_seg("alignA", 'chr1', 200, 100, "CCCNNNNNNN")
        actual_paired_alignment = connor.PairedAlignment(alignment1l, alignment1r)

        self.assertIs(alignment1l, actual_paired_alignment.left_alignment)
        self.assertIs(alignment1r, actual_paired_alignment.right_alignment)

    def test_eq(self):
        left = align_seg("alignA", 'chr1', 100, 200, "AAANNNNNNN")
        right = align_seg("alignA", 'chr1', 200, 100, "CCCNNNNNNN")
        other = align_seg("alignB", 'chr1', 100, 200, "AAANNNNNNN")

        base = connor.PairedAlignment(left, right)
        self.assertEquals(base, connor.PairedAlignment(left, right))
        self.assertNotEquals(base, connor.PairedAlignment(other, right))
        self.assertNotEquals(base, connor.PairedAlignment(left, other))

    def test_hash(self):
        left_A = align_seg("alignA", 'chr1', 100, 200, "AAANNNNNNN")
        right_A = align_seg("alignA", 'chr1', 200, 100, "CCCNNNNNNN")
        left_B = align_seg("alignA", 'chr1', 100, 200, "AAANNNNNNN")
        right_B = align_seg("alignA", 'chr1', 200, 100, "CCCNNNNNNN")

        actual_set = set()
        base = connor.PairedAlignment(left_A, right_A)
        actual_set.add(base)
        self.assertEquals(1, len(actual_set))

        actual_set.add(base)
        self.assertEquals(1, len(actual_set))

        actual_set.add(connor.PairedAlignment(left_A, right_A))
        self.assertEquals(1, len(actual_set))

        actual_set.add(connor.PairedAlignment(left_B, right_B))
        self.assertEquals(1, len(actual_set))


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
        align_A0 = align_seg("alignA", 'chr1', 10, 100)
        align_A1 = align_seg("alignA", 'chr1', 100, 10)
        align_B0 = align_seg("alignB", 'chr1', 10, 100)
        align_B1 = align_seg("alignB", 'chr1', 100, 10)
        alignments = [align_A0, align_B0, align_A1, align_B1]
        coord_read_name_manifest = {('chr1', 10, 100): set(['alignA', 'alignB'])}

        actual_families = [family for family in connor._build_coordinate_families(alignments, coord_read_name_manifest)]

        pair_A = connor.PairedAlignment(align_A0, align_A1)
        pair_B = connor.PairedAlignment(align_B0, align_B1)
        expected_families = [set([pair_A, pair_B])]
        self.assertEquals(expected_families, actual_families)

    def test_build_coordinate_families_threeFamilies(self):
        align_A0 = align_seg("alignA", 'chr1', 10, 100)
        align_A1 = align_seg("alignA", 'chr1', 100, 10)
        align_B0 = align_seg("alignB", 'chr1', 10, 100)
        align_B1 = align_seg("alignB", 'chr1', 100, 10)
        align_C0 = align_seg("alignC", 'chr1', 20, 200)
        align_C1 = align_seg("alignC", 'chr1', 200, 20)
        align_D0 = align_seg("alignD", 'chr1', 30, 300)
        align_D1 = align_seg("alignD", 'chr1', 300, 30)
        alignments = [align_A0, align_B0, align_C0, align_A1, align_B1,
                      align_D0, align_D1, align_C1]
        coord_read_name_manifest = {('chr1', 10, 100): set(['alignA', 'alignB']),
                                ('chr1', 20, 200): set(['alignC']),
                                ('chr1', 30, 300): set(['alignD'])}

        actual_families = [family for family in connor._build_coordinate_families(alignments, coord_read_name_manifest)]

        pair_A = connor.PairedAlignment(align_A0, align_A1)
        pair_B = connor.PairedAlignment(align_B0, align_B1)
        pair_C = connor.PairedAlignment(align_C0, align_C1)
        pair_D = connor.PairedAlignment(align_D0, align_D1)

        expected_families = [set([pair_A, pair_B]),
                             set([pair_D]),
                             set([pair_C])]
        self.assertEquals(expected_families, actual_families)

    def test_build_consensus_pair(self):
        align_A0 = align_seg("alignA", 'chr1', 10, 100)
        align_A1 = align_seg("alignA", 'chr1', 100, 10)
        align_B0 = align_seg("alignB", 'chr1', 10, 100)
        align_B1 = align_seg("alignB", 'chr1', 100, 10)
        alignments = set([connor.PairedAlignment(align_A0, align_A1),
                          connor.PairedAlignment(align_B0, align_B1)])

        actual_pair = connor._build_consensus_pair(alignments)

        expected_pair = connor.PairedAlignment(align_A0, align_A1)
        self.assertEquals(expected_pair, actual_pair)

    def test_build_tag_families(self):
        alignment1l = align_seg("alignA", 'chr1', 100, 200, "AAANNNNNNN")
        alignment1r = align_seg("alignA", 'chr1', 200, 100, "CCCNNNNNNN")
        alignment2l = align_seg("alignB", 'chr1', 100, 200, "GGGNNNNNNN")
        alignment2r = align_seg("alignB", 'chr1', 200, 100, "TTTNNNNNNN")
        alignment3l = align_seg("alignC", 'chr1', 100, 200, "AAANNNNNNN")
        alignment3r = align_seg("alignC", 'chr1', 200, 100, "CCCNNNNNNN")

        align_pair1 = connor.PairedAlignment(alignment1l, alignment1r)
        align_pair2 = connor.PairedAlignment(alignment2l, alignment2r)
        align_pair3 = connor.PairedAlignment(alignment3l, alignment3r)

        expected_tag_family_1 = frozenset([align_pair1, align_pair3])
        expected_tag_family_2 = frozenset([align_pair2])
        expected_tag_families = set([expected_tag_family_1,
                                     expected_tag_family_2])

        actual_tag_families = connor._build_tag_families(set([align_pair1,
                                                              align_pair2,
                                                              align_pair3]))

        self.assertEquals(expected_tag_families, actual_tag_families)


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


class ConnorFunctionalTestCase(unittest.TestCase):
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
            connor.main(input_bam, output_bam)
            alignments = _pysam_alignments_from_bam(output_bam)
            self.assertEquals(4, len(alignments))
            self.assertEquals(("readNameA1", 100),
                              (alignments[0].qname,
                               alignments[0].reference_start + 1))
            self.assertEquals(("readNameA1", 300),
                              (alignments[1].qname,
                               alignments[1].reference_start + 1))
            self.assertEquals(("readNameB1", 200),
                              (alignments[2].qname,
                               alignments[2].reference_start + 1))
            self.assertEquals(("readNameB1", 400),
                              (alignments[3].qname,
                               alignments[3].reference_start + 1))

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
            connor.main(input_bam, output_bam)
            alignments = _pysam_alignments_from_bam(output_bam)
            self.assertEquals(4, len(alignments))
            self.assertEquals(("readNameA1", 100),
                              (alignments[0].qname,
                               alignments[0].reference_start + 1))
            self.assertEquals(("readNameA1", 300),
                              (alignments[1].qname,
                               alignments[1].reference_start + 1))
            self.assertEquals(("readNameB1", 100),
                              (alignments[2].qname,
                               alignments[2].reference_start + 1))
            self.assertEquals(("readNameB1", 500),
                              (alignments[3].qname,
                               alignments[3].reference_start + 1))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
