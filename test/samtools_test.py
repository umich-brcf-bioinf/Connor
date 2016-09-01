#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments,deprecated-method
from __future__ import print_function, absolute_import, division
from argparse import Namespace
from copy import deepcopy
import os
import pysam
from testfixtures.tempdirectory import TempDirectory

import connor
from connor.samtools import BamFlag
from connor.samtools import BamTag
from connor.samtools import ConnorAlign
from connor.samtools import filter_alignments
import connor.samtools as samtools
from test.utils_test import MicroMock
import test.utils_test as utils_test


class MockAlignWriter(object):
    def __init__(self):
        self._write_calls = []
        self.bam_file_path = "foo.bam"
        self._close_was_called = False

    def write(self, family, connor_align):
        self._write_calls.append((family, connor_align))

    #pylint: disable=unused-argument
    def close(self, log=None):
        self._close_was_called = True

def mock_align(**kwargs):
    a = pysam.AlignedSegment()
    a.query_name = "align1"
    a.flag = 99
    a.reference_id = 0
    a.reference_start = 32
    a.mapping_quality = 20
    a.cigar = ((0,3), (2,1), (0,1))
    a.next_reference_id = 0
    a.next_reference_start=199
    a.template_length=167
    #query_qualities must be set after query_sequence
    a.query_sequence=kwargs.pop('query_sequence', 'AGCTTAG')
    for (key, value) in kwargs.items():
        setattr(a, key, value)
    return a

def _create_file(path, filename, contents):
    filename = os.path.join(path, filename)
    with open(filename, 'wt') as new_file:
        new_file.write(contents)
        new_file.flush()
    return filename

def create_bam(path, filename, sam_contents, index=True):
    sam_filename = _create_file(path, filename, sam_contents)
    bam_filename = sam_filename.replace(".sam", ".bam")
    pysam_bam_from_sam(sam_filename, bam_filename, index)
    return bam_filename

def pysam_bam_from_sam(sam_filename, bam_filename, index=True):
    infile = samtools.alignment_file(sam_filename, "r")
    outfile = samtools.alignment_file(bam_filename, "wb", template=infile)
    for s in infile:
        outfile.write(s)
    infile.close()
    outfile.close()
    if index:
        samtools.index(bam_filename)

def pysam_alignments_from_bam(bam_filename):
    infile = samtools.alignment_file(bam_filename, "rb")
    aligned_segments = [s for s in infile]
    infile.close()
    return aligned_segments


class ConnorAlignTest(utils_test.BaseConnorTestCase):
    @staticmethod
    def byte_array_to_string(sequence):
        if isinstance(sequence, str):
            return sequence
        else:
            return str(sequence.decode("utf-8"))

    def test_eq(self):
        pysam_align = mock_align(query_name="align1")
        base =  ConnorAlign(pysam_align)
        self.assertEqual(base, base)
        self.assertEqual(base, ConnorAlign(pysam_align))
        self.assertEqual(base, ConnorAlign(mock_align(query_name = "align1")))
        different_pysam_align = ConnorAlign(mock_align(query_name = "align2"))
        self.assertNotEqual(base, different_pysam_align)
        different_filter = ConnorAlign(pysam_align)
        different_filter.filter_value = "foo; bar"
        self.assertNotEqual(base, different_filter)

    def test_gettersPassthroughToPysamAlignSegment(self):
        pysam_align = mock_align(query_name="queryname_1",
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
        pysam_align = mock_align(query_name="queryname_1",
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
        pysam_align = mock_align(query_name="queryname_1",
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
        pysam_align = mock_align(reference_start=100, next_reference_start=200)
        self.assertEqual('left', ConnorAlign(pysam_align).orientation)

    def test_orientation_right(self):
        pysam_align = mock_align(reference_start=200, next_reference_start=100)
        self.assertEqual('right', ConnorAlign(pysam_align).orientation)

    def test_orientation_samePositiveIsLeft(self):
        pysam_align = mock_align(flag=1,
                                 reference_start=100,
                                 next_reference_start=100)
        self.assertEqual('left', ConnorAlign(pysam_align).orientation)

    def test_orientation_sameNegativeIsRight(self):
        pysam_align = mock_align(flag=17,
                                 reference_start=100,
                                 next_reference_start=100)
        self.assertEqual('right', ConnorAlign(pysam_align).orientation)


class SamtoolsTest(utils_test.BaseConnorTestCase):
    def test_build_writer(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
@PG|ID:bwa|VN:0.5.5
@PG|ID:GATK|PN:foo|VN:1.0.3471
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>
'''.replace("|", "\t")

        with TempDirectory() as tmp_dir:
            input_bam = create_bam(tmp_dir.path,
                                   'input.sam',
                                   sam_contents)
            annotated_output_bam = os.path.join(tmp_dir.path, 'annotated.bam')
            tags = []
            args=Namespace(original_command_line=['command-line'],
                           simplify_pg_header=False)
            actual_writer = samtools.build_writer(input_bam,
                                                  annotated_output_bam,
                                                  tags,
                                                  args)
            actual_writer.close()

            actual_output = samtools.alignment_file(annotated_output_bam, 'rb',)
            expected_header = {'HD': {'GO': 'none',
                                      'SO': 'coordinate',
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
            self.assertEqual(expected_header, actual_output.header)

    def test_build_annotated_aligns_writer_nullIfNotSpecified(self):
        actual_writer = samtools.build_writer(input_bam='foo',
                                              output_bam='',
                                              tags=[],
                                              args=Namespace())
        self.assertEqual(samtools.AlignWriter.NULL, actual_writer)

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
            bam = create_bam(tmp_dir.path,
                             "input.sam",
                             sam_contents,
                              index=False)
            samtools.sort_and_index_bam(bam)
            alignments = samtools.alignment_file(bam, "rb").fetch()
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
                bam = create_bam(os.path.join(tmp_dir.path, "tmp"),
                                 "input.sam",
                                 sam_contents,
                                 index=False)
                bam_filename = os.path.basename(bam)

                samtools.sort_and_index_bam(os.path.join("tmp", bam_filename))

                alignments = samtools.alignment_file(bam, "rb").fetch()
                aligns = [(a.query_name,
                           a.reference_start + 1) for a in alignments]
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

    def test_filter_alignments_passthorughIncludedAligns(self):
        align1 = mock_align(query_name="align1")
        base = [align1]
        excluded_writer = MockAlignWriter()

        aligns = [align for align in filter_alignments(base,
                                                       excluded_writer)]

        self.assertEqual([ConnorAlign(align1)],aligns)
        self.assertEqual(0, len(excluded_writer._write_calls))

    def test_filter_alignments_skipsWriteIfNoExcludedWriter(self):
        flag = 99
        align1 = mock_align(query_name="align1", flag=flag)
        align2 = mock_align(query_name="align2", flag=flag^BamFlag.PROPER_PAIR)
        align3 = mock_align(query_name="align3", flag=flag)
        base = [align1, align2, align3]

        query_names = [x.query_name for x in filter_alignments(base)]

        self.assertEqual(["align1", "align3"], query_names)


    def test_filter_alignments_excludesUnpairedAligns(self):
        flag = 99
        align1 = mock_align(query_name="align1", flag=flag)
        align2 = mock_align(query_name="align2", flag=flag^BamFlag.PROPER_PAIR)
        align3 = mock_align(query_name="align3", flag=flag)
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
        align1 = mock_align(query_name="align1", flag=flag)
        align2 = mock_align(query_name="align2", flag=flag | BamFlag.SECONDARY)
        align3 = mock_align(query_name="align3", flag=flag)
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

    def test_filter_alignments_excludesQCFails(self):
        flag = 99
        align1 = mock_align(query_name="align1", flag=flag)
        align2 = mock_align(query_name="align2", flag=flag | BamFlag.QCFAIL)
        align3 = mock_align(query_name="align3", flag=flag)
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
        align1 = mock_align(query_name="align1", mapping_quality=1)
        align2 = mock_align(query_name="align2", mapping_quality=0)
        align3 = mock_align(query_name="align3", mapping_quality=1)
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
        align1 = mock_align(query_name="align1", cigarstring="6M")
        align2 = mock_align(query_name="align2", cigarstring="*")
        align3 = mock_align(query_name="align3", cigarstring="6M")
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


class AlignWriterTest(utils_test.BaseConnorTestCase):
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
            align1 = ConnorAlign(mock_align(query_name="align1"))
            align2 = ConnorAlign(mock_align(query_name="align2"))
            align3 = ConnorAlign(mock_align(query_name="align3"))
            family = None
            writer = samtools.AlignWriter(header, bam_path)

            writer.write(family, align1)
            writer.write(family, align2)
            writer.write(family, align3)
            writer.close()

            bamfile = samtools.alignment_file(bam_path, 'rb')
            actual_query_names = [align.query_name for align in bamfile.fetch()]
            bamfile.close()

        self.assertEqual(['align1', 'align2', 'align3'], actual_query_names)

    def test_write_addsAlignTags(self):
        with TempDirectory() as tmp_dir:
            bam_path = os.path.join(tmp_dir.path, 'destination.bam')
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}] }
            align1 = ConnorAlign(mock_align(query_name='align1'))
            align2 = ConnorAlign(mock_align(query_name='align2'))
            align3 = ConnorAlign(mock_align(query_name='align3'))

            tag1 = BamTag('X1','Z', 'desc',
                          get_value=lambda family, align: family)
            tag2 = BamTag('X2','Z', 'desc',
                          get_value=lambda family, align: align.query_name)

            writer = samtools.AlignWriter(header, bam_path, [tag1, tag2])

            writer.write('familyA', align1)
            writer.write('familyB', align2)
            writer.write('familyC', align3)
            writer.close()

            bamfile = samtools.alignment_file(bam_path, 'rb')
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
        self.assertEqual("X2:Z:align1", align_tags[('align1', 'X2')])
        self.assertEqual("X2:Z:align2", align_tags[('align2', 'X2')])
        self.assertEqual("X2:Z:align3", align_tags[('align3', 'X2')])

    def test_write_skipsTagsWhenValueIsNone(self):
        with TempDirectory() as tmp_dir:
            bam_path = os.path.join(tmp_dir.path, 'destination.bam')
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}] }
            align1 = ConnorAlign(mock_align(query_name='align1'))
            align2 = ConnorAlign(mock_align(query_name='align2'))
            align3 = ConnorAlign(mock_align(query_name='align3'))

            get_value = lambda family, align: 'Yes' if align.query_name == 'align2' else None
            tag1 = BamTag('X1','Z', 'desc', get_value=get_value)

            writer = samtools.AlignWriter(header, bam_path, [tag1])

            writer.write('familyA', align1)
            writer.write('familyB', align2)
            writer.write('familyC', align3)
            writer.close()

            bamfile = samtools.alignment_file(bam_path, 'rb')
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
            align1 = ConnorAlign(mock_align(query_name='align1'))
            align1.set_tag('X1', 'No', 'Z')

            tag1 = BamTag('X1','Z', 'desc',
                          get_value = lambda family, align: None)

            writer = samtools.AlignWriter(header, bam_path, [tag1])

            writer.write('familyA', align1)
            writer.close()

            bamfile = samtools.alignment_file(bam_path, 'rb')
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
            tag1 = BamTag('X1','Z', 'annotates family', get_value=None)
            tag2 = BamTag('X2','Z', 'annotates alignment', get_value=None)
            writer = samtools.AlignWriter(header, bam_path, [tag2, tag1])
            writer.close()

            bamfile = samtools.alignment_file(bam_path, 'rb')
            actual_header = dict(bamfile.header)
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
        samtools.AlignWriter.NULL.write('family', 'connor_align')
        samtools.AlignWriter.NULL.close()
        self.assertEqual(1, 1)

    def test_close_sortsAndIndexes(self):
        with TempDirectory() as tmp_dir:
            bam_path = os.path.join(tmp_dir.path, 'destination.bam')
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}] }
            align1 = ConnorAlign(mock_align(query_name='align1',
                                            reference_start=100))
            align2 = ConnorAlign(mock_align(query_name='align2',
                                            reference_start=200))
            align3 = ConnorAlign(mock_align(query_name='align3',
                                            reference_start=300))

            tag1 = BamTag('X1','Z', 'desc',
                          get_value=lambda family, align: family)
            tag2 = BamTag('X2','Z', 'desc',
                          get_value=lambda family, align: align.query_name)

            writer = samtools.AlignWriter(header, bam_path, [tag1, tag2])

            writer.write('familyC', align3)
            writer.write('familyA', align1)
            writer.write('familyB', align2)
            writer.close()

            bamfile = samtools.alignment_file(bam_path, 'rb')
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
            align1 = ConnorAlign(mock_align(query_name='align1',
                                            reference_start=100))

            writer = samtools.AlignWriter(header, bam_path, [])

            writer.write('familyA', align1)
            writer.close(log=self.mock_logger)
        info_log_lines = self.mock_logger._log_calls['INFO']
        self.assertEqual(1, len(info_log_lines))
        self.assertRegexpMatches(info_log_lines[0], 'destination.bam')

class BamTagTest(utils_test.BaseConnorTestCase):
    def test_init_setsHeaderComment(self):
        tag = BamTag('foo', 'Z', 'foo description', lambda fam, align: None)
        self.assertEqual('connor\tBAM tag\tfoo: foo description',
                         tag.header_comment)

    def test_set_tag(self):
        get_value = lambda family, align: family + ":" + align.query_name
        tag = BamTag('X9', 'Z', 'foo description', get_value)
        connor_align = ConnorAlign(mock_align())

        tag.set_tag('family1', connor_align)

        self.assertEqual([('X9', 'family1:align1')], connor_align.get_tags())

    def test_lt_sortsByNameThenDescription(self):
        base = BamTag('X2', 'i', 'Desc B', None)
        self.assertEqual(False, base.__lt__(base))
        self.assertEqual(False, base.__lt__(BamTag('X2','i', 'Desc B', None)))

        self.assertEqual(True, base.__lt__(BamTag('X2','i', 'Desc C', None)))
        self.assertEqual(True, base.__lt__(BamTag('X3','i', 'Desc B', None)))

        self.assertEqual(False, base.__lt__(BamTag('X1','i', 'Desc B', None)))
        self.assertEqual(False, base.__lt__(BamTag('X2','i', 'Desc A', None)))


class LoggingWriterTest(utils_test.BaseConnorTestCase):
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
            writer.write(family, align)
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


    def test_close_whenAllPlaced(self):
        base_writer = MockAlignWriter()
        writer = samtools.LoggingWriter(base_writer, self.mock_logger)
        fam1 = MicroMock(umi_sequence=4, filter_value=None)
        alignA = MicroMock(filter_value=None)
        family_aligns = [(fam1, alignA), (fam1, alignA)]

        for family, align in family_aligns:
            writer.write(family, align)
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


    def test_write_passThroughToBaseWriter(self):
        base_writer = MockAlignWriter()
        writer = samtools.LoggingWriter(base_writer, self.mock_logger)
        fam1 = MicroMock(umi_sequence = 1, filter_value=None)
        al1A = MicroMock(filter_value=None)
        al1B = MicroMock(filter_value = 'foo')
        family_aligns = [(fam1, al1A), (fam1, al1B)]

        for family, align in family_aligns:
            writer.write(family, align)

        self.assertEqual([(fam1, al1A), (fam1, al1B)],
                         base_writer._write_calls)

    def test_write_UnplacedAlignWritesFamilyNone(self):
        base_writer = MockAlignWriter()
        writer = samtools.LoggingWriter(base_writer, self.mock_logger)
        fam1 = samtools.LoggingWriter.UNPLACED_FAMILY
        al1A = MicroMock(filter_value = 'foo')

        writer.write(fam1, al1A)

        self.assertEqual([(None, al1A)], base_writer._write_calls)

    def test_close_passThroughToBaseWriter(self):
        base_writer = MockAlignWriter()
        writer = samtools.LoggingWriter(base_writer, self.mock_logger)

        writer.close()

        self.assertEqual(True, base_writer._close_was_called)
