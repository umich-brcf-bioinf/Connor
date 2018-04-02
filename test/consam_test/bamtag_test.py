#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments,deprecated-method
from __future__ import print_function, absolute_import, division

from connor.consam import bamtag
from connor.consam.bamtag import BamTag
from connor.samtools import ConnorAlign
from test.utils_test import MicroMock
from test.utils_test import BaseConnorTestCase


class BamTagPackageTest(BaseConnorTestCase):
    @staticmethod
    def get_tag(tags, name):
        for tag in tags:
            if tag._tag_name == name:
                return tag
        return None

    def test_build_bam_tags(self):
        actual_tags = bamtag.build_bam_tags()
        self.assertEqual(7, len(actual_tags))

    def test_build_bam_tags_x0_filter(self):
        tag = BamTagPackageTest.get_tag(bamtag.build_bam_tags(), 'X0')
        self.assertEqual('X0', tag._tag_name)
        self.assertEqual('Z', tag._tag_type)
        self.assertRegexpMatches(tag._description, 'filter')

        self.assertEqual(None, tag._get_value(None, None, None))

        family = MicroMock(filter_value=None)
        connor_align = MicroMock(filter_value=None)
        self.assertEqual(None, tag._get_value(family, None, connor_align))

        family = MicroMock(filter_value='foo')
        connor_align = MicroMock(filter_value='bar')
        self.assertEqual('foo', tag._get_value(family, None, None))
        self.assertEqual('bar', tag._get_value(None, None, connor_align))
        self.assertEqual('foo;bar', tag._get_value(family, None, connor_align))

    def test_build_bam_tags_x1_positions(self):
        tag = BamTagPackageTest.get_tag(bamtag.build_bam_tags(), 'X1')
        self.assertEqual('Z', tag._tag_type)
        self.assertRegexpMatches(tag._description,
                                 'leftmost~rightmost matched pair positions')

        align = self.mock_align()
        pair = MicroMock(positions=lambda x:'100~150')
        tag.set_tag(None, pair, align)
        self.assertEqual([('X1', '100~150')], align.get_tags())

        align = self.mock_align()
        tag.set_tag(None, None, align)
        self.assertEqual([], align.get_tags())


    def test_build_bam_tags_x2_cigars(self):
        tag = BamTagPackageTest.get_tag(bamtag.build_bam_tags(), 'X2')
        self.assertEqual('Z', tag._tag_type)
        self.assertRegexpMatches(tag._description,
                                 'L~R CIGARs')

        align = self.mock_align()
        pair = MicroMock(cigars=lambda x:'1S2M4S~8S16M32S')
        tag.set_tag(None, pair, align)
        self.assertEqual([('X2', '1S2M4S~8S16M32S')], align.get_tags())

        align = self.mock_align()
        tag.set_tag(None, None, align)
        self.assertEqual([], align.get_tags())

    def test_build_bam_tags_x3_unique_identifier(self):
        tag = BamTagPackageTest.get_tag(bamtag.build_bam_tags(), 'X3')
        self.assertEqual('i', tag._tag_type)
        self.assertRegexpMatches(tag._description, 'unique identifier')

        family = MicroMock(umi_sequence=42)
        align = self.mock_align()
        tag.set_tag(family, None, align)
        self.assertEqual([('X3', 42)], align.get_tags())

        align = self.mock_align()
        tag.set_tag(None, None, align)
        self.assertEqual([], align.get_tags())

    def test_build_bam_tags_x4_umt_barcodes(self):
        tag = BamTagPackageTest.get_tag(bamtag.build_bam_tags(), 'X4')
        self.assertEqual('Z', tag._tag_type)
        self.assertRegexpMatches(tag._description, 'UMT barcodes')

        family = MicroMock(umt=lambda *args:'AAA~CCC')

        align = self.mock_align()
        tag.set_tag(family, None, align)
        self.assertEqual([('X4', 'AAA~CCC')], align.get_tags())

        align = self.mock_align()
        tag.set_tag(None, None, align)
        self.assertEqual([], align.get_tags())

    def test_build_bam_tags_x5_family_size(self):
        tag = BamTagPackageTest.get_tag(bamtag.build_bam_tags(), 'X5')
        self.assertEqual('i', tag._tag_type)
        self.assertRegexpMatches(tag._description, 'family size')

        family = MicroMock(included_pair_count=42)
        align = self.mock_align()
        tag.set_tag(family, None, align)
        self.assertEqual([('X5', 42)], align.get_tags())

        align = self.mock_align()
        tag.set_tag(None, None, align)
        self.assertEqual([], align.get_tags())

    def test_build_bam_tags_x6_consensus_template(self):
        tag = BamTagPackageTest.get_tag(bamtag.build_bam_tags(), 'X6')
        self.assertEqual('i', tag._tag_type)
        self.assertRegexpMatches(tag._description,
                                 'template for the consensus alignment')

        align = self.mock_align()
        family = MicroMock(is_consensus_template=lambda x: 1)
        tag.set_tag(family, None, align)
        self.assertEqual([('X6', 1)], align.get_tags())

        align = self.mock_align()
        family = MicroMock(is_consensus_template=lambda x: None)
        tag.set_tag(family, None, align)
        self.assertEqual([], align.get_tags())

        align = self.mock_align()
        tag.set_tag(None, None, align)
        self.assertEqual([], align.get_tags())

class BamTagTest(BaseConnorTestCase):
    def test_init_setsHeaderComment(self):
        tag = BamTag('foo', 'Z', 'foo description', lambda fam, align: None)
        self.assertEqual('connor\tBAM tag\tfoo: foo description',
                         tag.header_comment)

    def test_set_tag(self):
        def get_value (family, pair, align):
            return family + ':' + pair + ':' + align.query_name
        tag = BamTag('X9', 'Z', 'foo description', get_value)
        connor_align = ConnorAlign(self.mock_align())

        tag.set_tag('family1', 'pair1', connor_align)

        self.assertEqual([('X9', 'family1:pair1:align1')],
                         connor_align.get_tags())

    def test_set_tag_NoneReplacedWIthNullObject(self):
        def get_value(family, pair, align):
            return ':'.join([type(family).__name__,
                             type(pair).__name__,
                             type(align).__name__])
        tag = BamTag('X9', 'Z', 'foo description', get_value)
        connor_align = ConnorAlign(self.mock_align(query_name='baz'))

        tag.set_tag(None, None, connor_align)

        self.assertEqual([('X9', '_NullObject:_NullObject:ConnorAlign')], connor_align.get_tags())


    def test_lt_sortsByNameThenDescription(self):
        base = BamTag('X2', 'i', 'Desc B', None)
        self.assertEqual(False, base.__lt__(base))
        self.assertEqual(False, base.__lt__(BamTag('X2','i', 'Desc B', None)))

        self.assertEqual(True, base.__lt__(BamTag('X2','i', 'Desc C', None)))
        self.assertEqual(True, base.__lt__(BamTag('X3','i', 'Desc B', None)))

        self.assertEqual(False, base.__lt__(BamTag('X1','i', 'Desc B', None)))
        self.assertEqual(False, base.__lt__(BamTag('X2','i', 'Desc A', None)))
