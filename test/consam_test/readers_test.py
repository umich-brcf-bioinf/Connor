#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments
#pylint: disable=deprecated-method
from __future__ import print_function, absolute_import, division

from test.consam_test.writers_test import MockAlignWriter
from test.utils_test import BaseConnorTestCase

from connor.consam.alignments import ConnorAlign
from connor.consam.bamflag import BamFlag
import connor.consam.readers as readers

class ReadersTest(BaseConnorTestCase):
    def test_filter_alignments_passthorughIncludedAligns(self):
        align1 = self.mock_align(query_name="align1")
        base = [align1]
        excluded_writer = MockAlignWriter()

        aligns = [align for align in readers._filter_alignments(base,
                                                       excluded_writer)]

        self.assertEqual([ConnorAlign(align1)],aligns)
        self.assertEqual(0, len(excluded_writer._write_calls))

    def test_filter_alignments_skipsWriteIfNoExcludedWriter(self):
        flag = 99
        align1 = self.mock_align(query_name="align1", flag=flag)
        align2 = self.mock_align(query_name="align2", flag=flag^BamFlag.PROPER_PAIR)
        align3 = self.mock_align(query_name="align3", flag=flag)
        base = [align1, align2, align3]

        query_names = [x.query_name for x in readers._filter_alignments(base)]

        self.assertEqual(["align1", "align3"], query_names)


    def test_filter_alignments_excludesUnpairedAligns(self):
        flag = 99
        align1 = self.mock_align(query_name="align1", flag=flag)
        align2 = self.mock_align(query_name="align2", flag=flag^BamFlag.PROPER_PAIR)
        align3 = self.mock_align(query_name="align3", flag=flag)
        base = [align1, align2, align3]
        excluded_writer = MockAlignWriter()

        query_names = [x.query_name for x in readers._filter_alignments(base,
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

        names = [x.query_name for x in readers._filter_alignments(base,
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

        names = [x.query_name for x in readers._filter_alignments(base,
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

        names = [x.query_name for x in readers._filter_alignments(base,
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

        names = [x.query_name for x in readers._filter_alignments(base,
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

        names = [x.query_name for x in readers._filter_alignments(base,
                                                         excluded_writer)]

        self.assertEqual(["align1", "align3"], names)
        self.assertEqual(1, len(excluded_writer._write_calls))
        (family, connor_align) = excluded_writer._write_calls[0]
        self.assertEqual(None, family)
        self.assertEqual(align2.query_name, connor_align.query_name)
        self.assertEqual('cigar unavailable', connor_align.filter_value)


    def test_build_coordinate_pairs_singlePair(self):
        align1L = ConnorAlign(self.mock_align(query_name='1',
                                         reference_start=100,
                                         next_reference_start=200))
        align1R = ConnorAlign(self.mock_align(query_name='1',
                                         reference_start=200,
                                         next_reference_start=100))
        aligns = [align1L, align1R]

        umt_length = 6
        actual_pairs = [p for p in readers._build_coordinate_pairs(umt_length, aligns, None)]

        self.assertEqual(1, len(actual_pairs))
        actual_pair = actual_pairs[0]
        self.assertEqual(align1L, actual_pair.left)
        self.assertEqual(align1R, actual_pair.right)

    def test_build_coordinate_pairs_oneCoordinate(self):
        align1L = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=200))
        align2L = ConnorAlign(self.mock_align(query_name = '2',
                                         reference_start=100,
                                         next_reference_start=200))
        align1R = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=200,
                                         next_reference_start=100))
        align2R = ConnorAlign(self.mock_align(query_name = '2',
                                         reference_start=200,
                                         next_reference_start=100))
        aligns = [align1L, align2L, align1R, align2R]

        umt_length = 6
        actual_pairs = [f for f in readers._build_coordinate_pairs(umt_length, aligns, None)]

        self.assertEqual(2, len(actual_pairs))
        pairs = dict([(pair.query_name, pair) for pair in actual_pairs])
        self.assertEqual(align1L, pairs['1'].left)
        self.assertEqual(align1R, pairs['1'].right)
        self.assertEqual(align2L, pairs['2'].left)
        self.assertEqual(align2R, pairs['2'].right)

    def test_build_coordinate_pairs_twoCoordinatesSameRight(self):
        align1L = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=200))
        align2L = ConnorAlign(self.mock_align(query_name = '2',
                                         reference_start=100,
                                         next_reference_start=200))
        align3L = ConnorAlign(self.mock_align(query_name = '3',
                                         reference_start=125,
                                         next_reference_start=200))
        align1R = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=200,
                                         next_reference_start=100))
        align2R = ConnorAlign(self.mock_align(query_name = '2',
                                         reference_start=200,
                                         next_reference_start=100))
        align3R = ConnorAlign(self.mock_align(query_name = '3',
                                         reference_start=200,
                                         next_reference_start=125))
        aligns = [align1L, align2L, align3L, align1R, align2R, align3R]

        umt_length = 6
        actual_pairs = [f for f in readers._build_coordinate_pairs(umt_length, aligns, None)]

        self.assertEqual(3, len(actual_pairs))
        actual_pair_names = set([pair.query_name for pair in actual_pairs])
        self.assertEqual(set(['1', '2', '3']), actual_pair_names)

    def test_build_coordinate_pairs_identicalCoordinates(self):
        align1L = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=100))
        align1R = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=100))
        aligns = [align1L, align1R]
        writer = MockAlignWriter()
        umt_length = 6
        actual_pairs = [f for f in readers._build_coordinate_pairs(umt_length, aligns, writer)]

        self.assertEqual(1, len(actual_pairs))
        actual_pair_names = set([pair.query_name for pair in actual_pairs])
        self.assertEqual(set(['1']), actual_pair_names)


    def test_build_coordinate_pairs_orphanedRightIsSafe(self):
        align1L = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=200))
        align1R = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=200,
                                         next_reference_start=100))
        align2R = ConnorAlign(self.mock_align(query_name = '2',
                                         reference_start=200,
                                         next_reference_start=100))
        aligns = [align1L, align1R, align2R]
        writer = MockAlignWriter()
        umt_length = 6
        actual_pairs = [f for f in readers._build_coordinate_pairs(umt_length, aligns, writer)]

        self.assertEqual(1, len(actual_pairs))
        actual_pair_names = set([pair.query_name for pair in actual_pairs])
        self.assertEqual(set(['1']), actual_pair_names)

    def test_build_coordinate_pairs_orphanedRightWrittenToExcluded(self):
        align1L = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=200))
        align2R = ConnorAlign(self.mock_align(query_name = '2',
                                         reference_start=150,
                                         next_reference_start=100))
        align1R = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=200,
                                         next_reference_start=100))
        aligns = [align1L, align1R, align2R]
        writer = MockAlignWriter()

        umt_length = 6
        for _ in readers._build_coordinate_pairs(umt_length, aligns, writer):
            pass

        self.assertEqual(1, len(writer._write_calls))
        (actual_family, actual_align) = writer._write_calls[0]
        self.assertEqual(None, actual_family)
        self.assertEqual(align2R, actual_align)
        self.assertEqual('read mate was missing or excluded',
                         actual_align.filter_value)

    def test_build_coordinate_pairs_whenExhasutedRemaindersWrittenToExcluded(self):
        align1L = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=200))
        align1R = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=200,
                                         next_reference_start=100))
        align3L = ConnorAlign(self.mock_align(query_name = '3',
                                         reference_start=150,
                                         next_reference_start=250))
        aligns = [align1L, align3L, align1R]
        writer = MockAlignWriter()
        umt_length = 6
        for _ in readers._build_coordinate_pairs(umt_length, aligns, writer):
            pass

        self.assertEqual(1, len(writer._write_calls))
        (actual_family, actual_align) = writer._write_calls[0]
        self.assertEqual(None, actual_family)
        self.assertEqual(align3L, actual_align)
        self.assertEqual('read mate was missing or excluded',
                         actual_align.filter_value)

    def test_build_coordinate_pairs_lookForPassedPops(self):
        align1L = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=100,
                                         next_reference_start=200))
        align2L = ConnorAlign(self.mock_align(query_name = '2',
                                         reference_start=100,
                                         next_reference_start=200))
        align3L = ConnorAlign(self.mock_align(query_name = '3',
                                         reference_start=300,
                                         next_reference_start=400))
        align1R = ConnorAlign(self.mock_align(query_name = '1',
                                         reference_start=200,
                                         next_reference_start=100))
        align2R = ConnorAlign(self.mock_align(query_name = '2',
                                         reference_start=200,
                                         next_reference_start=100))

        class BombError(ValueError):
            pass

        class BombAlign(object):
            @property
            def orientation(self):
                raise BombError("Boom.")
        align4L = BombAlign()

        aligns = [align1L, align2L, align1R, align2R, align3L, align4L]
        actual_pairs = []
        umt_length = 6
        try:
            for pair in readers._build_coordinate_pairs(umt_length, aligns, None):
                actual_pairs.append(pair)
            self.fail("expected BombError")
        except BombError:
            pass

        self.assertEqual(2, len(actual_pairs))
        actual_pair_names = set([pair.query_name for pair in actual_pairs])
        self.assertEqual(set(['1', '2']), actual_pair_names)

    def test_progress_logger_passesThroughItems(self):
        base_gen = [0,1,2,3,4,5,6,7,8,9]
        total_rows = len(base_gen)
        actual_items = []
        for item in readers._progress_logger(base_gen,
                                         total_rows,
                                         self.mock_logger):
            actual_items.append(item)
        self.assertEqual(actual_items, base_gen)

    def test_progress_logger_logsProgress(self):
        base_gen = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
        total_rows = len(base_gen)
        gen = readers._progress_logger(base_gen,
                                     total_rows,
                                     self.mock_logger)
        next(gen)
        self.assertEqual("0% (1/20) alignments processed",
                         self.mock_logger._log_calls['INFO'][0])
        next(gen)
        self.assertEqual("10% (2/20) alignments processed",
                         self.mock_logger._log_calls['INFO'][1])
        next(gen)
        next(gen)
        self.assertEqual("20% (4/20) alignments processed",
                         self.mock_logger._log_calls['INFO'][2])
        next(gen)
        next(gen)
        self.assertEqual("30% (6/20) alignments processed",
                         self.mock_logger._log_calls['INFO'][3])
        for _ in range(7, 21):
            next(gen)

        self.assertRaises(StopIteration,
                          next,
                          gen)

        self.assertEqual("100% (20/20) alignments processed",
                         self.mock_logger._log_calls['INFO'][-1])

    def test_progress_logger_logsMem(self):
        base_gen = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
        total_rows = len(base_gen)
        def supplemental_log():
            self.mock_logger.debug("foo")
            self.mock_logger.debug("bar")
        gen = readers._progress_logger(base_gen,
                                     total_rows,
                                     self.mock_logger,
                                     supplemental_log)
        next(gen)
        next(gen)
        self.assertRegexpMatches(self.mock_logger._log_calls['DEBUG'][0],
                                 "foo")
        self.assertRegexpMatches(self.mock_logger._log_calls['DEBUG'][1],
                                 "bar")
        next(gen)
        self.assertRegexpMatches(self.mock_logger._log_calls['DEBUG'][2],
                                 "foo")
        self.assertRegexpMatches(self.mock_logger._log_calls['DEBUG'][3],
                                 "bar")
        next(gen)
        self.assertRegexpMatches(self.mock_logger._log_calls['DEBUG'][4],
                                 "foo")
        self.assertRegexpMatches(self.mock_logger._log_calls['DEBUG'][5],
                                 "bar")
