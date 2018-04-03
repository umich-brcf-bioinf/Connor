#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments,deprecated-method
from __future__ import print_function, absolute_import, division
from argparse import Namespace
from collections import OrderedDict
from copy import deepcopy
import os

from testfixtures.tempdirectory import TempDirectory

from connor import __version__
from connor.consam.alignments import ConnorAlign
import connor.consam.bamtag as bamtag
import connor.consam.pysamwrapper as pysamwrapper
import connor.consam.writers as writers

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


class WritersTest(BaseConnorTestCase):
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
            actual_writer = writers.build_writer(input_bam,
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
                                       'VN':__version__,
                                       'CL':'command-line'
                                       },
                                      ]}
            actual_header = pysamwrapper.get_header_dict(actual_output)
            if isinstance(actual_header, OrderedDict):
                actual_header = dict(actual_header)
            self.assertEqual(expected_header, actual_header)

    def test_build_annotated_aligns_writer_nullIfNotSpecified(self):
        actual_writer = writers.build_writer(input_bam='foo',
                                              output_bam='',
                                              tags=[],
                                              args=Namespace())
        self.assertEqual(writers.AlignWriter.NULL, actual_writer)


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
            writer = writers.AlignWriter(header, bam_path)
            writer.close()
        self.assertEqual(bam_path, writer._bam_path)

    def test_init_defaultToNoTags(self):
        with TempDirectory() as tmp_dir:
            bam_path = os.path.join(tmp_dir.path, "destination.bam")
            header = { 'HD': {'VN': '1.0'},
                      'SQ': [{'LN': 1575, 'SN': 'chr1'},
                             {'LN': 1584, 'SN': 'chr2'}] }
            writer = writers.AlignWriter(header, bam_path)
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
            writer = writers.AlignWriter(header, bam_path)

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

            writer = writers.AlignWriter(header, bam_path, [tag1, tag2, tag3])

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

            writer = writers.AlignWriter(header, bam_path, [tag1])

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

            writer = writers.AlignWriter(header, bam_path, [tag1])

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
            writer = writers.AlignWriter(header, bam_path, [tag2, tag1])
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
        writers.AlignWriter.NULL.write('family', None, 'connor_align')
        writers.AlignWriter.NULL.close()
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

            writer = writers.AlignWriter(header, bam_path, [tag1, tag2])

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

            writer = writers.AlignWriter(header, bam_path, [])

            writer.write('familyA', None, align1)
            writer.close(log=self.mock_logger)
        info_log_lines = self.mock_logger._log_calls['INFO']
        self.assertEqual(1, len(info_log_lines))
        self.assertRegexpMatches(info_log_lines[0], 'destination.bam')


class LoggingWriterTest(BaseConnorTestCase):
    def test_close_logsFilterStats(self):
        base_writer = MockAlignWriter()
        writer = writers.LoggingWriter(base_writer, self.mock_logger)
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
        writer = writers.LoggingWriter(base_writer, self.mock_logger)
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
        writer = writers.LoggingWriter(base_writer, self.mock_logger)
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
        writer = writers.LoggingWriter(base_writer, self.mock_logger)
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
        writer = writers.LoggingWriter(base_writer, self.mock_logger)

        included_count, total_count, filter_counts = writer._family_stats

        self.assertEqual(0, included_count)
        self.assertEqual(0, total_count)
        self.assertEqual({}, filter_counts)

    def test_write_passThroughToBaseWriter(self):
        base_writer = MockAlignWriter()
        writer = writers.LoggingWriter(base_writer, self.mock_logger)
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
        writer = writers.LoggingWriter(base_writer, self.mock_logger)
        fam1 = writers.LoggingWriter.UNPLACED_FAMILY
        al1A = MicroMock(filter_value = 'foo')

        writer.write(fam1, None, al1A)

        self.assertEqual([(None, al1A)], base_writer._write_calls)

    def test_close_passThroughToBaseWriter(self):
        base_writer = MockAlignWriter()
        writer = writers.LoggingWriter(base_writer, self.mock_logger)

        writer.close()

        self.assertEqual(True, base_writer._close_was_called)
