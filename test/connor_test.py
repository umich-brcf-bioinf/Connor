#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring
from __future__ import print_function, absolute_import
from collections import namedtuple
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


class TestConnor(unittest.TestCase):

    def test_build_alignment_family_dict(self):
        Align = namedtuple('Align', 'name key')
        align1 = Align(name='align1', key=3)
        align2 = Align(name='align2', key=3)
        align3 = Align(name='align3', key=4)

        actual_dict = connor._build_alignment_family_dict([align1,
                                                           align2,
                                                           align3])

        expected_dict = {3: set(['align1', 'align2']),
                         4: set(['align3'])}
        self.assertEquals(expected_dict, actual_dict)

    def test_build_read_families_oneFamily(self):
        align_A0 = align_seg("alignA", 'chr1', 10, 100)
        align_A1 = align_seg("alignA", 'chr1', 100, 10)
        align_B0 = align_seg("alignB", 'chr1', 10, 100)
        align_B1 = align_seg("alignB", 'chr1', 100, 10)
        alignments = [align_A0, align_B0, align_A1, align_B1]
        coord_read_name_dict = {('chr1', 10, 100): set(['alignA', 'alignB'])}

        actual_families = [family for family in connor._build_read_families(alignments, coord_read_name_dict)]

        expected_families = [set([align_A0, align_A1, align_B0, align_B1])]
        self.assertEquals(expected_families, actual_families)

    def test_build_read_families_threeFamilies(self):
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
        coord_read_name_dict = {('chr1', 10, 100): set(['alignA', 'alignB']),
                                ('chr1', 20, 200): set(['alignC']),
                                ('chr1', 30, 300): set(['alignD'])}

        actual_families = [family for family in connor._build_read_families(alignments, coord_read_name_dict)]

        expected_families = [set([align_A0, align_A1, align_B0, align_B1]),
                             set([align_D0, align_D1]),
                             set([align_C0, align_C1])]
        self.assertEquals(expected_families, actual_families)

    def test_build_consensus_read(self):
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




if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
