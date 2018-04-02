#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments,deprecated-method
from __future__ import print_function, absolute_import, division
from argparse import Namespace
from collections import OrderedDict
from copy import deepcopy
import os


import pysam
from testfixtures.tempdirectory import TempDirectory

import connor
from connor.consam.bamflag import BamFlag
import connor.consam.pysamwrapper as pysamwrapper
import connor.consam.bamtag as bamtag
import connor.samtools as samtools
from connor.samtools import ConnorAlign
from connor.samtools import filter_alignments
from test.utils_test import MicroMock
from test.utils_test import BaseConnorTestCase


class MockAlignWriter(object):
    def __init__(self):
        self._write_calls = []
        self.bam_file_path = "foo.bam"
        self._close_was_called = False

    #pylint: disable=unused-argument
    def write(self, family, paired_align, connor_align):
        self._write_calls.append((family, connor_align))

    #pylint: disable=unused-argument
    def close(self, log=None):
        self._close_was_called = True

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
        actual_paired_alignment = samtools.PairedAlignment(left_align,
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
        paired_alignment = samtools.PairedAlignment(left, right, tag_length=1)
        self.assertEqual((101,251), paired_alignment.positions())
        self.assertEqual('101~251',
                         paired_alignment.positions('{left}~{right}'))

    def test_filter_value(self):
        left = ConnorAlign(self.mock_align(), filter_value=None)
        right = ConnorAlign(self.mock_align(), filter_value=None)
        paired_alignment = samtools.PairedAlignment(left, right, tag_length=1)
        self.assertEqual(None, paired_alignment.filter_value)

        left = ConnorAlign(self.mock_align(), filter_value='')
        right = ConnorAlign(self.mock_align(), filter_value='')
        paired_alignment = samtools.PairedAlignment(left, right, tag_length=1)
        self.assertEqual(None, paired_alignment.filter_value)

        left = ConnorAlign(self.mock_align(), filter_value='foo')
        right = ConnorAlign(self.mock_align(), filter_value=None)
        paired_alignment = samtools.PairedAlignment(left, right, tag_length=1)
        self.assertEqual(('foo', None), paired_alignment.filter_value)

        left = ConnorAlign(self.mock_align(), filter_value=None)
        right = ConnorAlign(self.mock_align(), filter_value='bar')
        paired_alignment = samtools.PairedAlignment(left, right, tag_length=1)
        self.assertEqual((None, 'bar'), paired_alignment.filter_value)

    def test_query_name(self):
        left = self.mock_align(query_name="alignA", reference_start=100)
        right = self.mock_align(query_name="alignA", reference_start=200)
        paired_alignment = samtools.PairedAlignment(left, right, tag_length=1)
        self.assertEqual("alignA", paired_alignment.query_name)

    def test_eq(self):
        left = self.mock_align(reference_start=100, next_reference_start=200)
        right = self.mock_align(reference_start=200, next_reference_start=100)
        other = self.mock_align(reference_start=0, next_reference_start=500)

        base = samtools.PairedAlignment(left, right)
        self.assertEquals(base, samtools.PairedAlignment(left, right))
        self.assertNotEquals(base, samtools.PairedAlignment(other, right))
        self.assertNotEquals(base, samtools.PairedAlignment(left, other))

    def test_hash(self):
        left_A = self.mock_align(query_name="alignA", reference_start=100)
        right_A = self.mock_align(query_name="alignA", reference_start=200)
        left_B = self.mock_align(query_name="alignA", reference_start=100)
        right_B = self.mock_align(query_name="alignA", reference_start=200)

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
        left_A = self.mock_align(query_sequence='AANN', query_qualities=[1,2,3,4])
        right_A = self.mock_align(query_sequence='NNCC', query_qualities=[5,6,7,8])
        paired_align = samtools.PairedAlignment(left_A, right_A, tag_length=2)

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


class SamtoolsTest(BaseConnorTestCase):
    def test_total_align_count(self):
        self.check_sysout_safe()
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameB1|147|chr10|400|0|5M|=|200|100|CCCCC|>>>>>
readNameA1|147|chr10|300|0|5M|=|100|100|AAAAA|>>>>>
readNameA1|99|chr10|100|0|5M|=|300|200|AAAAA|>>>>>
readNameB1|99|chr10|200|0|5M|=|400|200|CCCCC|>>>>>
readNameC1|12|chr10|400|0|*|=|200|100|CCCCC|>>>>>
readNameC1|12|chr10|400|0|*|=|200|100|CCCCC|>>>>>
readNameZ1|77|*|0|0|*|*|0|0|TTTTT|>>>>>
readNameZ1|141|*|0|0|*|*|0|0|GGGGG|>>>>>
'''.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam = self.create_bam(tmp_dir.path,
                                       'input.sam',
                                       sam_contents,
                                       index=False)
            pysamwrapper.sort_and_index_bam(input_bam)
            actual_count = samtools.total_align_count(input_bam)
            self.assertEqual(6, actual_count)

    def test_build_writer(self):
        sam_contents = \
'''@HD|VN:1.4|SO:coordinate
@SQ|SN:chr10|LN:135534747
@PG|ID:bwa|VN:0.5.5
@PG|ID:GATK|PN:foo|VN:1.0.3471
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>
'''.replace("|", "\t")

        with TempDirectory() as tmp_dir:
            input_bam = self.create_bam(tmp_dir.path,
                                        'input.sam',
                                        sam_contents)
            annotated_output_bam = os.path.join(tmp_dir.path, 'annotated.bam')
            tags = []
            args = Namespace(original_command_line=['command-line'],
                             simplify_pg_header=False)
            actual_writer = samtools.build_writer(input_bam,
                                                  annotated_output_bam,
                                                  tags,
                                                  args)
            actual_writer.close()

            actual_output = pysamwrapper.alignment_file(annotated_output_bam, 'rb',)
            expected_header = {'HD': {'SO': 'coordinate',
                                      'VN': '1.4'},
                               'SQ': [{'SN': 'chr10', 'LN': 135534747}],
                               'PG': [{'ID':'bwa', 'VN':'0.5.5'},
                                      {'ID':'GATK', 'PN':'foo', 'VN':'1.0.3471'},
                                      {'ID':'connor',
                                       'PN':'connor',
                                       'VN':connor.__version__,
                                       'CL':'command-line'
                                       },
                                      ]}
            actual_header = pysamwrapper.get_header_dict(actual_output)
            if isinstance(actual_header, OrderedDict):
                actual_header = dict(actual_header)
            self.assertEqual(expected_header, actual_header)

    def test_build_annotated_aligns_writer_nullIfNotSpecified(self):
        actual_writer = samtools.build_writer(input_bam='foo',
                                              output_bam='',
                                              tags=[],
                                              args=Namespace())
        self.assertEqual(samtools.AlignWriter.NULL, actual_writer)

    def test_filter_alignments_passthorughIncludedAligns(self):
        align1 = self.mock_align(query_name="align1")
        base = [align1]
        excluded_writer = MockAlignWriter()

        aligns = [align for align in filter_alignments(base,
                                                       excluded_writer)]

        self.assertEqual([ConnorAlign(align1)],aligns)
        self.assertEqual(0, len(excluded_writer._write_calls))

    def test_filter_alignments_skipsWriteIfNoExcludedWriter(self):
        flag = 99
        align1 = self.mock_align(query_name="align1", flag=flag)
        align2 = self.mock_align(query_name="align2", flag=flag^BamFlag.PROPER_PAIR)
        align3 = self.mock_align(query_name="align3", flag=flag)
        base = [align1, align2, align3]

        query_names = [x.query_name for x in filter_alignments(base)]

        self.assertEqual(["align1", "align3"], query_names)


    def test_filter_alignments_excludesUnpairedAligns(self):
        flag = 99
        align1 = self.mock_align(query_name="align1", flag=flag)
        align2 = self.mock_align(query_name="align2", flag=flag^BamFlag.PROPER_PAIR)
        align3 = self.mock_align(query_name="align3", flag=flag)
        base = [align1, align2, align3]
        excluded_writer = MockAlignWriter()

        query_names = [x.query_name for x in filter_alignments(base,
                                                               excluded_writer)]

        self.assertEqual(["align1", "align3"], query_names)
        self.assertEqual(1, len(excluded_writer._write_calls))
        (family, connor_align) = excluded_writer._write_calls[0]
        self.assertEqual(None, family)
        self.assertEqual(align2.query_name, connor_align.query_name)
        self.assertEqual('not in proper pair', connor_align.filter_value)

    def test_filter_alignments_excludesSecondaryAligns(self):
        flag = 99
        align1 = self.mock_align(query_name="align1", flag=flag)
        align2 = self.mock_align(query_name="align2", flag=flag | BamFlag.SECONDARY)
        align3 = self.mock_align(query_name="align3", flag=flag)
        base = [align1, align2, align3]
        excluded_writer = MockAlignWriter()

        names = [x.query_name for x in filter_alignments(base,
                                                         excluded_writer)]

        self.assertEqual(["align1", "align3"], names)
        self.assertEqual(1, len(excluded_writer._write_calls))
        (family, connor_align) = excluded_writer._write_calls[0]
        self.assertEqual(None, family)
        self.assertEqual(align2.query_name, connor_align.query_name)
        self.assertEqual('secondary alignment', connor_align.filter_value)

    def test_filter_alignments_excludesSupplementaryAligns(self):
        flag = 99
        align1 = self.mock_align(query_name="align1", flag=flag)
        align2 = self.mock_align(query_name="align2",
                            flag=flag | BamFlag.SUPPLEMENTARY)
        align3 = self.mock_align(query_name="align3", flag=flag)
        base = [align1, align2, align3]
        excluded_writer = MockAlignWriter()

        names = [x.query_name for x in filter_alignments(base,
                                                         excluded_writer)]

        self.assertEqual(["align1", "align3"], names)
        self.assertEqual(1, len(excluded_writer._write_calls))
        (family, connor_align) = excluded_writer._write_calls[0]
        self.assertEqual(None, family)
        self.assertEqual(align2.query_name, connor_align.query_name)
        self.assertEqual('supplementary alignment', connor_align.filter_value)


    def test_filter_alignments_excludesQCFails(self):
        flag = 99
        align1 = self.mock_align(query_name="align1", flag=flag)
        align2 = self.mock_align(query_name="align2", flag=flag | BamFlag.QCFAIL)
        align3 = self.mock_align(query_name="align3", flag=flag)
        base = [align1, align2, align3]
        excluded_writer = MockAlignWriter()

        names = [x.query_name for x in filter_alignments(base,
                                                         excluded_writer)]

        self.assertEqual(["align1", "align3"], names)
        self.assertEqual(1, len(excluded_writer._write_calls))
        (family, connor_align) = excluded_writer._write_calls[0]
        self.assertEqual(None, family)
        self.assertEqual(align2.query_name, connor_align.query_name)
        self.assertEqual('qc failed', connor_align.filter_value)

    def test_filter_alignments_excludesMapq0(self):
        align1 = self.mock_align(query_name="align1", mapping_quality=1)
        align2 = self.mock_align(query_name="align2", mapping_quality=0)
        align3 = self.mock_align(query_name="align3", mapping_quality=1)
        base = [align1, align2, align3]
        excluded_writer = MockAlignWriter()

        names = [x.query_name for x in filter_alignments(base,
                                                         excluded_writer)]

        self.assertEqual(["align1", "align3"], names)
        self.assertEqual(1, len(excluded_writer._write_calls))
        (family, connor_align) = excluded_writer._write_calls[0]
        self.assertEqual(None, family)
        self.assertEqual(align2.query_name, connor_align.query_name)
        self.assertEqual('mapping quality < 1', connor_align.filter_value)

    def test_filter_alignments_excludesCigarUnavailable(self):
        align1 = self.mock_align(query_name="align1", cigarstring="6M")
        align2 = self.mock_align(query_name="align2", cigarstring="*")
        align3 = self.mock_align(query_name="align3", cigarstring="6M")
        base = [align1, align2, align3]
        excluded_writer = MockAlignWriter()

        names = [x.query_name for x in filter_alignments(base,
                                                         excluded_writer)]

        self.assertEqual(["align1", "align3"], names)
        self.assertEqual(1, len(excluded_writer._write_calls))
        (family, connor_align) = excluded_writer._write_calls[0]
        self.assertEqual(None, family)
        self.assertEqual(align2.query_name, connor_align.query_name)
        self.assertEqual('cigar unavailable', connor_align.filter_value)


class AlignWriterTest(BaseConnorTestCase):
    @staticmethod
    def fix_pysam_inconsistent_tag_type(t_type):
        try:
            t_type = chr(t_type)
        except TypeError:
            pass
        return t_type

    def test_init_saves_file_path(self):
        with TempDirectory() as tmp_dir:
            bam_path = os.path.join(tmp_dir.path, "destination.bam")
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}] }
            writer = samtools.AlignWriter(header, bam_path)
            writer.close()
        self.assertEqual(bam_path, writer._bam_path)

    def test_init_defaultToNoTags(self):
        with TempDirectory() as tmp_dir:
            bam_path = os.path.join(tmp_dir.path, "destination.bam")
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}] }
            writer = samtools.AlignWriter(header, bam_path)
            writer.close()
        self.assertEqual([], writer._tags)

    def test_write(self):
        with TempDirectory() as tmp_dir:
            bam_path = os.path.join(tmp_dir.path, "destination.bam")
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}] }
            align1 = ConnorAlign(self.mock_align(query_name="align1"))
            align2 = ConnorAlign(self.mock_align(query_name="align2"))
            align3 = ConnorAlign(self.mock_align(query_name="align3"))
            family = None
            writer = samtools.AlignWriter(header, bam_path)

            writer.write(family, None, align1)
            writer.write(family, None, align2)
            writer.write(family, None, align3)
            writer.close()

            bamfile = pysamwrapper.alignment_file(bam_path, 'rb')
            actual_query_names = [align.query_name for align in bamfile.fetch()]
            bamfile.close()

        self.assertEqual(['align1', 'align2', 'align3'], actual_query_names)

    def test_write_addsAlignTags(self):
        with TempDirectory() as tmp_dir:
            bam_path = os.path.join(tmp_dir.path, 'destination.bam')
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}] }
            align1 = ConnorAlign(self.mock_align(query_name='align1'))
            align2 = ConnorAlign(self.mock_align(query_name='align2'))
            align3 = ConnorAlign(self.mock_align(query_name='align3'))

            tag1 = bamtag.BamTag('X1','Z', 'desc',
                          get_value=lambda family,pair,align: family)
            tag2 = bamtag.BamTag('X2','Z', 'desc',
                          get_value=lambda family,pair,align: pair)
            tag3 = bamtag.BamTag('X3','Z', 'desc',
                          get_value=lambda family,pair,align: align.query_name)

            writer = samtools.AlignWriter(header, bam_path, [tag1, tag2, tag3])

            writer.write('familyA', 'pair1', align1)
            writer.write('familyB', 'pair2', align2)
            writer.write('familyC', 'pair3', align3)
            writer.close()

            bamfile = pysamwrapper.alignment_file(bam_path, 'rb')
            actual_aligns = [a for a in bamfile.fetch()]
            bamfile.close()

        align_tags = {}
        for actual_align in actual_aligns:
            for t_name, t_val, t_type  in actual_align.get_tags(with_value_type=True):
                key = (actual_align.query_name, t_name)
                t_type = AlignWriterTest.fix_pysam_inconsistent_tag_type(t_type)
                align_tags[key] = "{}:{}:{}".format(t_name, t_type, t_val)

        self.assertEqual(3, len(actual_aligns))
        self.assertEqual("X1:Z:familyA", align_tags[('align1', 'X1')])
        self.assertEqual("X1:Z:familyB", align_tags[('align2', 'X1')])
        self.assertEqual("X1:Z:familyC", align_tags[('align3', 'X1')])
        self.assertEqual("X2:Z:pair1", align_tags[('align1', 'X2')])
        self.assertEqual("X2:Z:pair2", align_tags[('align2', 'X2')])
        self.assertEqual("X2:Z:pair3", align_tags[('align3', 'X2')])
        self.assertEqual("X3:Z:align1", align_tags[('align1', 'X3')])
        self.assertEqual("X3:Z:align2", align_tags[('align2', 'X3')])
        self.assertEqual("X3:Z:align3", align_tags[('align3', 'X3')])

    def test_write_skipsTagsWhenValueIsNone(self):
        with TempDirectory() as tmp_dir:
            bam_path = os.path.join(tmp_dir.path, 'destination.bam')
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}] }
            align1 = ConnorAlign(self.mock_align(query_name='align1'))
            align2 = ConnorAlign(self.mock_align(query_name='align2'))
            align3 = ConnorAlign(self.mock_align(query_name='align3'))

            get_value = lambda family, pair, align: 'Yes' if align.query_name == 'align2' else None
            tag1 = bamtag.BamTag('X1','Z', 'desc', get_value=get_value)

            writer = samtools.AlignWriter(header, bam_path, [tag1])

            writer.write('familyA', None, align1)
            writer.write('familyB', None, align2)
            writer.write('familyC', None, align3)
            writer.close()

            bamfile = pysamwrapper.alignment_file(bam_path, 'rb')
            actual_aligns = [a for a in bamfile.fetch()]
            bamfile.close()

        align_tags = {}
        for actual_align in actual_aligns:
            for t_name, t_val, t_type  in actual_align.get_tags(with_value_type=True):
                key = (actual_align.query_name, t_name)
                t_type = AlignWriterTest.fix_pysam_inconsistent_tag_type(t_type)
                align_tags[key] = "{}:{}:{}".format(t_name, t_type, t_val)

        self.assertEqual(3, len(actual_aligns))
        self.assertEqual(1, len(align_tags))
        self.assertEqual("X1:Z:Yes", align_tags[('align2', 'X1')])


    def test_write_removesTagsWhenValueIsNone(self):
        with TempDirectory() as tmp_dir:
            bam_path = os.path.join(tmp_dir.path, 'destination.bam')
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}] }
            align1 = ConnorAlign(self.mock_align(query_name='align1'))
            align1.set_tag('X1', 'No', 'Z')

            tag1 = bamtag.BamTag('X1','Z', 'desc',
                          get_value = lambda family, pair, align: None)

            writer = samtools.AlignWriter(header, bam_path, [tag1])

            writer.write('familyA', None, align1)
            writer.close()

            bamfile = pysamwrapper.alignment_file(bam_path, 'rb')
            actual_aligns = [a for a in bamfile.fetch()]
            bamfile.close()

        align_tags = {}
        for actual_align in actual_aligns:
            for t_name, t_val, t_type  in actual_align.get_tags(with_value_type=True):
                key = (actual_align.query_name, t_name)
                t_type = AlignWriterTest.fix_pysam_inconsistent_tag_type(t_type)
                align_tags[key] = "{}:{}:{}".format(t_name, t_type, t_val)

        self.assertEqual(1, len(actual_aligns))
        self.assertEqual(0, len(align_tags))


    def test_write_addsHeaderTags(self):
        with TempDirectory() as tmp_dir:
            bam_path = os.path.join(tmp_dir.path, 'destination.bam')
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}],
                      'CO': ['comment1', 'comment2']}
            tag1 = bamtag.BamTag('X1','Z', 'annotates family', get_value=None)
            tag2 = bamtag.BamTag('X2','Z', 'annotates alignment', get_value=None)
            writer = samtools.AlignWriter(header, bam_path, [tag2, tag1])
            writer.close()

            bamfile = pysamwrapper.alignment_file(bam_path, 'rb')
            actual_header = dict(pysamwrapper.get_header_dict(bamfile))
            bamfile.close()

        expected_header = deepcopy(header)
        expected_header.pop('CO')
        actual_comments = actual_header.pop('CO')
        expected_comments = ['comment1',
                             'comment2',
                             'connor\tBAM tag\tX1: annotates family',
                             'connor\tBAM tag\tX2: annotates alignment']
        self.assertEqual(expected_comments, actual_comments)

    def test_null_writer_methods(self):
        samtools.AlignWriter.NULL.write('family', None, 'connor_align')
        samtools.AlignWriter.NULL.close()
        self.assertEqual(1, 1)

    def test_close_sortsAndIndexes(self):
        with TempDirectory() as tmp_dir:
            bam_path = os.path.join(tmp_dir.path, 'destination.bam')
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}] }
            align1 = ConnorAlign(self.mock_align(query_name='align1',
                                            reference_start=100))
            align2 = ConnorAlign(self.mock_align(query_name='align2',
                                            reference_start=200))
            align3 = ConnorAlign(self.mock_align(query_name='align3',
                                            reference_start=300))

            tag1 = bamtag.BamTag('X1','Z', 'desc',
                          get_value=lambda family, pair, align: family)
            tag2 = bamtag.BamTag('X2','Z', 'desc',
                          get_value=lambda family, pair, align: align.query_name)

            writer = samtools.AlignWriter(header, bam_path, [tag1, tag2])

            writer.write('familyC', None, align3)
            writer.write('familyA', None, align1)
            writer.write('familyB', None, align2)
            writer.close()

            bamfile = pysamwrapper.alignment_file(bam_path, 'rb')
            actual_aligns = [a for a in bamfile.fetch()]
            bamfile.close()

            self.assertEqual(3, len(actual_aligns))
            self.assertEqual('align1', actual_aligns[0].query_name)
            self.assertEqual('align2', actual_aligns[1].query_name)
            self.assertEqual('align3', actual_aligns[2].query_name)

    def test_close_logs(self):
        with TempDirectory() as tmp_dir:
            bam_path = os.path.join(tmp_dir.path, 'destination.bam')
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}] }
            align1 = ConnorAlign(self.mock_align(query_name='align1',
                                            reference_start=100))

            writer = samtools.AlignWriter(header, bam_path, [])

            writer.write('familyA', None, align1)
            writer.close(log=self.mock_logger)
        info_log_lines = self.mock_logger._log_calls['INFO']
        self.assertEqual(1, len(info_log_lines))
        self.assertRegexpMatches(info_log_lines[0], 'destination.bam')


class LoggingWriterTest(BaseConnorTestCase):
    def test_close_logsFilterStats(self):
        base_writer = MockAlignWriter()
        writer = samtools.LoggingWriter(base_writer, self.mock_logger)
        fam1 = None
        al1A = MicroMock(filter_value='low mapping qual')
        al1B = MicroMock(filter_value='low mapping qual')
        fam2 = None
        al2A = MicroMock(filter_value='unpaired read')
        al2B = MicroMock(filter_value='unpaired read')
        fam3 = MicroMock(umi_sequence=3, filter_value=None)
        al3A = MicroMock(filter_value='minority CIGAR')
        al3B = MicroMock(filter_value=None)
        fam4 = MicroMock(umi_sequence=4, filter_value=None)
        al4A = MicroMock(filter_value=None)
        al4B = MicroMock(filter_value=None)
        fam5 = MicroMock(umi_sequence=5, filter_value='small family')
        al5A = MicroMock(filter_value=None)
        al5B = MicroMock(filter_value=None)
        family_aligns = [(fam1, al1A), (fam1, al1B),
                         (fam2, al2A), (fam2, al2B),
                         (fam3, al3A), (fam3, al3B),
                         (fam4, al4A), (fam4, al4B),
                         (fam5, al5A), (fam5, al5B)]

        for family, align in family_aligns:
            writer.write(family, None, align)
        writer.close()

        log_lines = self.mock_logger._log_calls['INFO']
        self.assertEqual('70.00% (7/10) alignments unplaced or discarded',
                         log_lines[0])
        self.assertEqual('families discarded: 33.33% (1/3) small family',
                         log_lines[1])
        self.assertEqual('30.00% (3/10) alignments included in 2 families',
                         log_lines[2])
        self.assertEqual('33.33% deduplication rate (1 - 2 families/3 included alignments)',
                         log_lines[3])
        self.assertEqual(4, len(log_lines))

        log_lines = self.mock_logger._log_calls['DEBUG']
        self.assertEqual('alignments unplaced: 20.00% (2/10) low mapping qual',
                         log_lines[0])
        self.assertEqual('alignments unplaced: 20.00% (2/10) unpaired read',
                         log_lines[1])
        self.assertEqual('alignments discarded: 20.00% (2/10) small family',
                         log_lines[2])
        self.assertEqual('alignments discarded: 10.00% (1/10) minority CIGAR',
                         log_lines[3])
        self.assertEqual(4, len(log_lines))

    def test_close_logsFilterStatsWarnsWhenNoAlignments(self):
        base_writer = MockAlignWriter()
        writer = samtools.LoggingWriter(base_writer, self.mock_logger)
        fam1 = None
        al1A = MicroMock(filter_value='low mapping qual')
        al1B = MicroMock(filter_value='low mapping qual')
        fam2 = None
        al2A = MicroMock(filter_value='unpaired read')
        al2B = MicroMock(filter_value='unpaired read')
        family_aligns = [(fam1, al1A), (fam1, al1B),
                         (fam2, al2A), (fam2, al2B)]

        for family, align in family_aligns:
            writer.write(family, None, align)
        writer.close()

        log_lines = self.mock_logger._log_calls['WARNING']
        self.assertEqual('No alignments passed filters. (Was input BAM downsampled?)',
                         log_lines[0])
        self.assertEqual(1, len(log_lines))


    def test_close_whenAllPlaced(self):
        base_writer = MockAlignWriter()
        writer = samtools.LoggingWriter(base_writer, self.mock_logger)
        fam1 = MicroMock(umi_sequence=4, filter_value=None)
        alignA = MicroMock(filter_value=None)
        family_aligns = [(fam1, alignA), (fam1, alignA)]

        for family, align in family_aligns:
            writer.write(family, None, align)
        writer.close()

        log_lines = self.mock_logger._log_calls['INFO']
        self.assertEqual('0.00% (0/2) alignments unplaced or discarded',
                         log_lines[0])
        self.assertEqual('100.00% (2/2) alignments included in 1 families',
                         log_lines[1])
        self.assertEqual('50.00% deduplication rate (1 - 1 families/2 included alignments)',
                         log_lines[2])
        self.assertEqual(3, len(log_lines))

        log_lines = self.mock_logger._log_calls['DEBUG']
        self.assertEqual(0, len(log_lines))

    def test_family_stats(self):
        base_writer = MockAlignWriter()
        writer = samtools.LoggingWriter(base_writer, self.mock_logger)
        fam1 = MicroMock(umi_sequence=4, filter_value=None)
        alignA = MicroMock(filter_value=None)
        family_aligns = [(fam1, alignA), (fam1, alignA)]

        for family, align in family_aligns:
            writer.write(family, None, align)

        included_count, total_count, filter_counts = writer._family_stats

        self.assertEqual(1, included_count)
        self.assertEqual(1, total_count)
        self.assertEqual({}, filter_counts)

    def test_family_stats_noFamilies(self):
        base_writer = MockAlignWriter()
        writer = samtools.LoggingWriter(base_writer, self.mock_logger)

        included_count, total_count, filter_counts = writer._family_stats

        self.assertEqual(0, included_count)
        self.assertEqual(0, total_count)
        self.assertEqual({}, filter_counts)

    def test_write_passThroughToBaseWriter(self):
        base_writer = MockAlignWriter()
        writer = samtools.LoggingWriter(base_writer, self.mock_logger)
        fam1 = MicroMock(umi_sequence = 1, filter_value=None)
        al1A = MicroMock(filter_value=None)
        al1B = MicroMock(filter_value = 'foo')
        family_aligns = [(fam1, al1A), (fam1, al1B)]

        for family, align in family_aligns:
            writer.write(family, None, align)

        self.assertEqual([(fam1, al1A), (fam1, al1B)],
                         base_writer._write_calls)

    def test_write_UnplacedAlignWritesFamilyNone(self):
        base_writer = MockAlignWriter()
        writer = samtools.LoggingWriter(base_writer, self.mock_logger)
        fam1 = samtools.LoggingWriter.UNPLACED_FAMILY
        al1A = MicroMock(filter_value = 'foo')

        writer.write(fam1, None, al1A)

        self.assertEqual([(None, al1A)], base_writer._write_calls)

    def test_close_passThroughToBaseWriter(self):
        base_writer = MockAlignWriter()
        writer = samtools.LoggingWriter(base_writer, self.mock_logger)

        writer.close()

        self.assertEqual(True, base_writer._close_was_called)
