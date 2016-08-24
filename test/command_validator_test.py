#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments
#pylint: disable=deprecated-method

from __future__ import print_function, absolute_import, division
from argparse import Namespace
import os
from testfixtures.tempdirectory import TempDirectory
import connor.command_validator as validator
import connor.utils as utils
import test.samtools_test as samtools_test
from test.utils_test import BaseConnorTestCase


class MockTask(object):
    def __init__(self, error_message=None):
        self.error_message = error_message
        self.args = None
        self.execute_called = False
        self.log = None

    def execute(self, args, log):
        self.execute_called = True
        if self.error_message:
            raise Exception(self.error_message)
        else:
            self.args = args
            self.log = log


class CommandValidatorTest(BaseConnorTestCase):
    def test_Validations(self):
        function_names = [f.__name__ for f in validator._VALIDATIONS]
        self.assertEqual(['_check_input_bam_exists',
                          '_check_input_bam_valid',
                          '_check_input_bam_indexed',
                          '_check_input_bam_not_deduped',
                          '_check_input_bam_not_empty',
                          '_check_input_bam_paired'],
                         function_names)

    def test_preflight_runsAllValidations(self):
        task1 = MockTask()
        task2 = MockTask()
        validator._VALIDATIONS = [task1.execute,
                                  task2.execute]
        args = Namespace()
        log = self.mock_logger
        validator.preflight(args, log)

        self.assertTrue(task1.execute_called)
        self.assertEqual(task1.args, args)
        self.assertEqual(task1.log, log)
        self.assertTrue(task2.execute_called)
        self.assertEqual(task2.args, args)
        self.assertEqual(task2.log, log)


    def test_input_bam_exists(self):
        with TempDirectory() as tmp_dir:
            tmp_dir.write('input.bam', b'foo')
            input_bam_path = os.path.join(tmp_dir.path, 'input.bam')
            args = Namespace(input_bam=input_bam_path)
            validator._check_input_bam_exists(args)
            self.ok()

    def test_input_bam_exists_raisesUsageError(self):
        with TempDirectory() as tmp_dir:
            input_bam_path = os.path.join(tmp_dir.path, 'input.bam')
            args = Namespace(input_bam=input_bam_path)
            self.assertRaisesRegexp(utils.UsageError,
                                    r'\[.*input.bam\] does not exist',
                                    validator._check_input_bam_exists,
                                    args)

    def test_check_input_bam_valid(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>'''.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=False)
            args = Namespace(input_bam=input_bam_path)
            validator._check_input_bam_valid(args)
            self.ok()

    def test_check_input_bam_valid_raisesUsageError(self):
        with TempDirectory() as tmp_dir:
            tmp_dir.write('input.bam', b'foo')
            input_bam_path = os.path.join(tmp_dir.path, 'input.bam')
            args = Namespace(input_bam=input_bam_path)
            self.assertRaisesRegexp(utils.UsageError,
                                    r'\[.*input.bam\] not a valid BAM',
                                    validator._check_input_bam_valid,
                                    args)

    def test_check_input_bam_indexed(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>'''.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path)
            validator._check_input_bam_indexed(args)
            self.ok()


    def test_check_input_bam_indexed_raisesUsageError(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>'''.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=False)
            args = Namespace(input_bam=input_bam_path)
            self.assertRaisesRegexp(utils.UsageError,
                                    r'\[.*input.bam\] is not indexed',
                                    validator._check_input_bam_indexed,
                                    args)

    def test_check_input_bam_not_connor_generated(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
@PG|ID:foo|PN:bwa
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>'''.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=False)
            validator._check_input_bam_not_deduped(args)
            self.ok()


    def test_check_input_bam_not_connor_generated_raisesUsageError(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
@PG|ID:foo|PN:connor
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>'''.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=False)
            regex = (r'\[.*input.bam\] has already been processed with Connor'
                     r'.*Are you sure.*force')
            self.assertRaisesRegexp(utils.UsageError,
                                    regex,
                                    validator._check_input_bam_not_deduped,
                                    args)

    def test_check_input_bam_not_connor_generated_noPgHeader(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>'''.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=False)
            validator._check_input_bam_not_deduped(args)
            self.ok()

    def test_check_input_bam_not_connor_generated_noPnHeader(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
@PG|ID:bwa|VN:1.3
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>'''.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=False)
            validator._check_input_bam_not_deduped(args)
            self.ok()

    def test_check_input_bam_not_connor_generated_warnIfForce(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
@PG|ID:foo|PN:connor
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>'''.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=True)
            validator._check_input_bam_not_deduped(args,
                                                   log=self.mock_logger)
            warnings = self.mock_logger._log_calls['WARNING']
            self.assertEqual(1, len(warnings))
            regex = (r'\[.*input.bam\] has already been processed with Connor'
                     r'.*forcing')
            self.assertRegexpMatches(warnings[0], regex)

    def test_check_input_bam_not_empty_raiseUsageError(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
'''.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=False)
            self.assertRaisesRegexp(utils.UsageError,
                                    r'\[.*input.bam\] is empty',
                                    validator._check_input_bam_not_empty,
                                    args)

    def test_check_input_bam_not_empty(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|5M|=|300|200|AAAAA|>>>>>'''.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=False)
            validator._check_input_bam_not_empty(args)
            self.ok()

    def test_check_input_bam_paired_raisesUsageError(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|{flag}|chr10|100|20|5M|=|300|200|AAAAA|>>>>>'''
        sam_contents = sam_contents.format(flag='16').replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=False)
            regex = r'\[.*input.bam\] does not appear to contain paired reads'
            self.assertRaisesRegexp(utils.UsageError,
                                    regex,
                                    validator._check_input_bam_paired,
                                    args)

    def test_check_input_bam_paired_ok(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|{unpaired_flag}|chr10|100|20|5M|=|300|200|AAAAA|>>>>>
readNameA1|{paired_flag}|chr10|100|20|5M|=|300|200|AAAAA|>>>>>'''
        sam_contents = sam_contents.format(unpaired_flag='16', paired_flag='99')
        sam_contents = sam_contents.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=False)
            validator._check_input_bam_paired(args)
            self.ok()

    def test_check_input_bam_paired_forcingOk(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|{flag}|chr10|100|20|5M|=|300|200|AAAAA|>>>>>'''
        sam_contents = sam_contents.format(flag='16').replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=True)
            validator._check_input_bam_paired(args, self.mock_logger)
        warnings = self.mock_logger._log_calls['WARNING']
        self.assertEqual(1, len(warnings))
        regex = (r'\[.*input.bam\] does not appear to contain paired '
                 r'reads.*forcing')
        self.assertRegexpMatches(warnings[0], regex)


    def test_check_input_bam_barcoded_raisesUsageError_ok(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA1|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA2|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA2|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA3|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA3|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA4|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA4|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA5|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA5|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
'''
        sam_contents = sam_contents.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=False)
            validator._check_input_bam_barcoded(args)
            self.ok()

    def test_check_input_bam_barcoded_ok_at_threshold(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|8M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA1|147|chr10|100|20|8M|=|300|200|AAAAANNN|>>>>>>>>
readNameA2|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA2|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA3|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA3|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA4|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA4|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA5|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA5|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
'''
        sam_contents = sam_contents.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=False)
            validator._check_input_bam_barcoded(args)
            self.ok()

    def test_check_input_bam_barcoded_leftUnbarcodedRaisesUsageError(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|8M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA1|147|chr10|100|20|8M|=|300|200|AAAAANNN|>>>>>>>>
readNameA2|99|chr10|100|20|8M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA2|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA3|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA3|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA4|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA4|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA5|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA5|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
'''
        sam_contents = sam_contents.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=False)
            regex = r'\[.*input.bam\] reads do not appear to have barcodes'
            self.assertRaisesRegexp(utils.UsageError,
                                    regex,
                                    validator._check_input_bam_barcoded,
                                    args)

    def test_check_input_bam_barcoded_rightUnbarcodedRaisesUsageError(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|8M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA1|147|chr10|100|20|8M|=|300|200|AAAAANNN|>>>>>>>>
readNameA2|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA2|147|chr10|100|20|8M|=|300|200|AAAAANNN|>>>>>>>>
readNameA3|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA3|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA4|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA4|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA5|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA5|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
'''
        sam_contents = sam_contents.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=False)
            regex = r'\[.*input.bam\] reads do not appear to have barcodes'
            self.assertRaisesRegexp(utils.UsageError,
                                    regex,
                                    validator._check_input_bam_barcoded,
                                    args)

    def test_check_input_bam_barcoded_unbarcodedForcedOk(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameA1|99|chr10|100|20|8M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA1|147|chr10|100|20|8M|=|300|200|AAAAANNN|>>>>>>>>
readNameA2|99|chr10|100|20|8M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA2|147|chr10|100|20|8M|=|300|200|AAAAANNN|>>>>>>>>
readNameA3|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA3|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA4|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA4|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
readNameA5|99|chr10|100|20|3S5M|=|300|200|NNNAAAAA|>>>>>>>>
readNameA5|147|chr10|100|20|5M3S|=|300|200|AAAAANNN|>>>>>>>>
'''
        sam_contents = sam_contents.replace("|", "\t")
        with TempDirectory() as tmp_dir:
            input_bam_path = samtools_test.create_bam(tmp_dir.path,
                                                      "input.sam",
                                                      sam_contents,
                                                      index=True)
            args = Namespace(input_bam=input_bam_path, force=True)
            validator._check_input_bam_barcoded(args, self.mock_logger)
        warnings = self.mock_logger._log_calls['WARNING']
        self.assertEqual(1, len(warnings))
        regex = r'\[.*input.bam\] reads do not appear to have barcodes.*forcing'
        self.assertRegexpMatches(warnings[0], regex)
