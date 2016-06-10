#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring
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
                                      'next_reference_start'))

def align_seg(a, b, c, d):
    return MockPysamAlignedSegment(query_name=a,
                                   reference_name=b,
                                   reference_start=c,
                                   next_reference_start=d)

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

        expected_families = [set([align_A0, align_A1, align_B0, align_B1])]
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

        expected_families = [set([align_A0, align_A1, align_B0, align_B1]),
                             set([align_D0, align_D1]),
                             set([align_C0, align_C1])]
        self.assertEquals(expected_families, actual_families)

    def test_build_consensus_pair(self):
        align_A0 = align_seg("alignA", 'chr1', 10, 100)
        align_A1 = align_seg("alignA", 'chr1', 100, 10)
        align_B0 = align_seg("alignB", 'chr1', 10, 100)
        align_B1 = align_seg("alignB", 'chr1', 100, 10)
        alignments = set([align_A0, align_B0, align_B1, align_A1])

        actual_pair = connor._build_consensus_pair(alignments)

        expected_pair = set([align_B0, align_B1])
        self.assertEquals(expected_pair, set(actual_pair))


class TestLightweightAlignment(unittest.TestCase):
    def test_lightweight_alignment_forwardRead(self):
        alignedSegment = MockPysamAlignedSegment(query_name="align1",
                                                 reference_name='chr1',
                                                 reference_start=10,
                                                 next_reference_start=100)

        actual_lwa = connor.LightweightAlignment(alignedSegment)

        self.assertEquals("align1", actual_lwa.name)
        self.assertEquals(('chr1', 10, 100), actual_lwa.key)

    def test_lightweight_alignment_reverseRead(self):
        alignedSegment = MockPysamAlignedSegment(query_name="align1",
                                                 reference_name='chr1',
                                                 reference_start=100,
                                                 next_reference_start=10)

        actual_lwa = connor.LightweightAlignment(alignedSegment)

        self.assertEquals("align1", actual_lwa.name)
        self.assertEquals(('chr1', 10, 100), actual_lwa.key)

    def test_lightweight_alignment_weirdRead(self):
        alignedSegment = MockPysamAlignedSegment(query_name="align1",
                                                 reference_name='chr1',
                                                 reference_start=100,
                                                 next_reference_start=100)

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
            self.assertEquals(("readNameA2", 100),
                              (alignments[0].qname,
                               alignments[0].reference_start + 1))
            self.assertEquals(("readNameA2", 300),
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
