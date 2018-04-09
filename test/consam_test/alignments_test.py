#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments,deprecated-method
from __future__ import print_function, absolute_import, division

from test.utils_test import MicroMock
from test.utils_test import BaseConnorTestCase

from connor.consam.alignments import ConnorAlign
from connor.consam.alignments import PairedAlignment


class ConnorAlignTest(BaseConnorTestCase):
    def test_eq(self):
        pysam_align = self.mock_align(query_name="align1")
        base =  ConnorAlign(pysam_align)
        self.assertEqual(base, base)
        self.assertEqual(base, ConnorAlign(pysam_align))
        self.assertEqual(base, ConnorAlign(self.mock_align(query_name = "align1")))
        different_pysam_align = ConnorAlign(self.mock_align(query_name = "align2"))
        self.assertNotEqual(base, different_pysam_align)
        different_filter = ConnorAlign(pysam_align)
        different_filter.filter_value = "foo; bar"
        self.assertNotEqual(base, different_filter)

    def test_hash(self):
        pysam_align_A_42 = self.mock_align(query_name="A", reference_start=42)
        pysam_align_B_42 = self.mock_align(query_name="B", reference_start=42)
        pysam_align_A_43 = self.mock_align(query_name="A", reference_start=43)

        base = ConnorAlign(pysam_align_A_42, filter_value="f1")
        same = ConnorAlign(pysam_align_A_42, filter_value="f1")
        different_query_name = ConnorAlign(pysam_align_B_42, filter_value="f1")
        different_start = ConnorAlign(pysam_align_A_43, filter_value="f1")
        different_filter = ConnorAlign(pysam_align_A_42, filter_value="f2")
        self.assertEqual(base.__hash__(), same.__hash__())
        self.assertNotEqual(base.__hash__(), different_query_name.__hash__())
        self.assertNotEqual(base.__hash__(), different_start.__hash__())
        self.assertNotEqual(base.__hash__(), different_filter.__hash__())

    def test_gettersPassthroughToPysamAlignSegment(self):
        pysam_align = self.mock_align(query_name="queryname_1",
                            flag=99,
                            reference_id=3,
                            reference_start=142,
                            mapping_quality=20,
                            cigarstring="8M",
                            next_reference_id=4,
                            next_reference_start=242,
                            template_length=100,
                            query_sequence="ACGTACGT",
                            query_qualities=[20]*8,
                            )
        pysam_align.set_tag('X1', 'foo')
        connor_align = ConnorAlign(pysam_align)

        self.assertEqual('queryname_1', connor_align.query_name)
        self.assertEqual(99, connor_align.flag)
        self.assertEqual(3, connor_align.reference_id)
        self.assertEqual(142, connor_align.reference_start)
        self.assertEqual(20, connor_align.mapping_quality)
        self.assertEqual('8M', connor_align.cigarstring)
        self.assertEqual(242, connor_align.next_reference_start)
        self.assertEqual(100, connor_align.template_length)
        self.assertEqual('ACGTACGT',
                         ConnorAlignTest.byte_array_to_string(connor_align.query_sequence))
        self.assertEqual([20] * 8, connor_align.query_qualities)
        self.assertEqual(150, connor_align.reference_end)
        self.assertEqual('foo', connor_align.get_tag('X1'))
        self.assertEqual([('X1', 'foo')], connor_align.get_tags())


    def test_settersPassthroughToPysamAlignSegment(self):
        pysam_align = self.mock_align(query_name="queryname_1",
                            flag=99,
                            reference_id=3,
                            reference_start=142,
                            mapping_quality=20,
                            cigarstring="8M",
                            next_reference_id=4,
                            next_reference_start=242,
                            template_length=100,
                            query_sequence="ACGTACGT",
                            query_qualities=[20]*8,
                            )
        connor_align = ConnorAlign(pysam_align)


        connor_align.query_name = 'queryname_11'
        connor_align.flag = 147
        connor_align.reference_id = 13
        connor_align.reference_start = 1142
        connor_align.mapping_quality = 120
        connor_align.cigarstring = "2S8M"
        connor_align.next_reference_id = 14
        connor_align.next_reference_start = 1242
        connor_align.template_length = 1100
        connor_align.query_sequence = "TTACGTACGT"
        connor_align.query_qualities = [20]*10
        connor_align.set_tag('X1', 'foo', 'Z')

        self.assertEqual('queryname_11', pysam_align.query_name)
        self.assertEqual(147, pysam_align.flag)
        self.assertEqual(13, pysam_align.reference_id)
        self.assertEqual(1142, pysam_align.reference_start)
        self.assertEqual(120, pysam_align.mapping_quality)
        self.assertEqual('2S8M', pysam_align.cigarstring)
        self.assertEqual(1242, pysam_align.next_reference_start)
        self.assertEqual(1100, pysam_align.template_length)
        self.assertEqual('TTACGTACGT',
                         ConnorAlignTest.byte_array_to_string(pysam_align.query_sequence))
        self.assertEqual([20] * 10, pysam_align.query_qualities)
        self.assertEqual(1150, pysam_align.reference_end)
        self.assertEqual(('foo', 'Z'),
                         pysam_align.get_tag('X1', with_value_type=True))

    def test_filter(self):
        pysam_align = self.mock_align(query_name="queryname_1",
                            flag=99,
                            reference_id=3,
                            reference_start=142,
                            mapping_quality=20,
                            cigarstring="8M",
                            next_reference_id=4,
                            next_reference_start=242,
                            template_length=100,
                            query_sequence="ACGTACGT",
                            query_qualities=[20]*8,
                            )
        connor_align = ConnorAlign(pysam_align)

        self.assertEqual(None, connor_align.filter_value)
        connor_align.filter_value = 'foo'
        self.assertEqual('foo', connor_align.filter_value)

    def test_orientation_left(self):
        pysam_align = self.mock_align(reference_start=100, next_reference_start=200)
        self.assertEqual('left', ConnorAlign(pysam_align).orientation)

    def test_orientation_right(self):
        pysam_align = self.mock_align(reference_start=200, next_reference_start=100)
        self.assertEqual('right', ConnorAlign(pysam_align).orientation)

    def test_orientation_sameIsNeither(self):
        pysam_align = self.mock_align(flag=129,
                                 reference_start=100,
                                 next_reference_start=100)
        self.assertEqual('neither', ConnorAlign(pysam_align).orientation)


class PairedAlignmentTest(BaseConnorTestCase):
    def test_init(self):
        left_align = self.mock_align(query_name="alignA",
                                query_sequence="AAATTT" "GGGG")
        right_align = self.mock_align(query_name="alignA",
                                 query_sequence="TTTT" "CCCGGG")
        tag_length = 6
        actual_paired_alignment = PairedAlignment(left_align,
                                                         right_align,
                                                         tag_length)

        self.assertIs(left_align, actual_paired_alignment.left)
        self.assertIs(right_align, actual_paired_alignment.right)
        left_umt = self.byte_array_to_string(actual_paired_alignment.umt[0])
        right_umt = self.byte_array_to_string(actual_paired_alignment.umt[1])
        self.assertEquals(("AAATTT", "CCCGGG"), (left_umt, right_umt))

    def test_init_valueErrorOnInconsistentQueryNames(self):
        left = self.mock_align(query_name="alignA")
        right = self.mock_align(query_name="alignB")
        self.assertRaisesRegexp(ValueError,
                                (r'Inconsistent query names '
                                 r'\(alignA != alignB\)'),
                                PairedAlignment,
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
        paired_alignment = PairedAlignment(left, right, tag_length=1)
        self.assertEqual(('1S2M4S', '16S32M64S'), paired_alignment.cigars())
        self.assertEqual('1S2M4S~16S32M64S',
                         paired_alignment.cigars('{left}~{right}'))


    def test_positions(self):
        left = MicroMock(query_name='A',
                         reference_start=100,
                         reference_end=150,
                         query_sequence='AAAAAA')
        right = MicroMock(query_name='A',
                          reference_start=200,
                          reference_end=250,
                          query_sequence='AAAAAA')
        paired_alignment = PairedAlignment(left, right, tag_length=1)
        self.assertEqual((101,251), paired_alignment.positions())
        self.assertEqual('101~251',
                         paired_alignment.positions('{left}~{right}'))

    def test_filter_value(self):
        left = ConnorAlign(self.mock_align(), filter_value=None)
        right = ConnorAlign(self.mock_align(), filter_value=None)
        paired_alignment = PairedAlignment(left, right, tag_length=1)
        self.assertEqual(None, paired_alignment.filter_value)

        left = ConnorAlign(self.mock_align(), filter_value='')
        right = ConnorAlign(self.mock_align(), filter_value='')
        paired_alignment = PairedAlignment(left, right, tag_length=1)
        self.assertEqual(None, paired_alignment.filter_value)

        left = ConnorAlign(self.mock_align(), filter_value='foo')
        right = ConnorAlign(self.mock_align(), filter_value=None)
        paired_alignment = PairedAlignment(left, right, tag_length=1)
        self.assertEqual(('foo', None), paired_alignment.filter_value)

        left = ConnorAlign(self.mock_align(), filter_value=None)
        right = ConnorAlign(self.mock_align(), filter_value='bar')
        paired_alignment = PairedAlignment(left, right, tag_length=1)
        self.assertEqual((None, 'bar'), paired_alignment.filter_value)

    def test_query_name(self):
        left = self.mock_align(query_name="alignA", reference_start=100)
        right = self.mock_align(query_name="alignA", reference_start=200)
        paired_alignment = PairedAlignment(left, right, tag_length=1)
        self.assertEqual("alignA", paired_alignment.query_name)

    def test_eq(self):
        umt_length = 6
        left = self.mock_align(reference_start=100, next_reference_start=200)
        right = self.mock_align(reference_start=200, next_reference_start=100)
        other = self.mock_align(reference_start=0, next_reference_start=500)

        base = PairedAlignment(left, right, umt_length)
        self.assertEquals(base, PairedAlignment(left, right, umt_length))
        self.assertNotEquals(base, PairedAlignment(other, right, umt_length))
        self.assertNotEquals(base, PairedAlignment(left, other, umt_length))
        self.assertNotEquals(base, PairedAlignment(left, right, 1))

    def test_hash(self):
        umt_length = 6
        left_A = self.mock_align(query_name="alignA", reference_start=100)
        right_A = self.mock_align(query_name="alignA", reference_start=200)
        left_B = self.mock_align(query_name="alignA", reference_start=100)
        right_B = self.mock_align(query_name="alignA", reference_start=200)

        actual_set = set()
        base = PairedAlignment(left_A, right_A, umt_length)
        actual_set.add(base)
        self.assertEquals(1, len(actual_set))

        actual_set.add(base)
        self.assertEquals(1, len(actual_set))

        actual_set.add(PairedAlignment(left_A, right_A, umt_length))
        self.assertEquals(1, len(actual_set))

        equivalent_pair = PairedAlignment(left_B, right_B, umt_length)
        actual_set.add(equivalent_pair)
        self.assertEquals(1, len(actual_set))

    def test_replace_umt(self):
        left_A = self.mock_align(query_sequence='AANN', query_qualities=[1,2,3,4])
        right_A = self.mock_align(query_sequence='NNCC', query_qualities=[5,6,7,8])
        paired_align = PairedAlignment(left_A, right_A, tag_length=2)

        paired_align.replace_umt(('GG','TT'))

        left = paired_align.left
        right = paired_align.right
        self.assertEquals('GGNN',
                          self.byte_array_to_string(left.query_sequence))
        self.assertEquals('NNTT',
                          self.byte_array_to_string(right.query_sequence))
        self.assertEquals([1,2,3,4],
                          left.query_qualities)
        self.assertEquals([5,6,7,8],
                          right.query_qualities)

    def test_replace_umt_errorIfInconsistentUmtLength(self):
        left_A = self.mock_align(query_sequence='AANN', query_qualities=[1,2,3,4])
        right_A = self.mock_align(query_sequence='NNCC', query_qualities=[5,6,7,8])
        paired_align = PairedAlignment(left_A, right_A, tag_length=2)

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
