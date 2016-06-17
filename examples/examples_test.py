from __future__ import print_function, absolute_import
import filecmp
import os
import unittest

from testfixtures import TempDirectory

import connor.connor as connor

INPUT_DIR=os.path.realpath(os.path.dirname(__file__))

class MockLogger(object):
    def __init__(self):
        self._log_calls = []

    def log(self, msg_format, *args):
        self._log_calls.append((msg_format, args))


class ExamplesFunctionalTest(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.original_logger = connor._log
        self.mock_logger = MockLogger()
        connor._log = self.mock_logger.log

    def tearDown(self):
        connor._log = self.original_logger
        unittest.TestCase.tearDown(self)
    
    def test_examples(self):
        with TempDirectory() as output_dir:

            input_bam = "PIK3CA-original.bam"
            output_bam = "PIK3CA-deduped.bam"
            output_bai = output_bam + ".bai"
            input_bam_filename = os.path.join(INPUT_DIR, input_bam)
            expect_bam_filename = os.path.join(INPUT_DIR, output_bam)
            expect_bai_filename = os.path.join(INPUT_DIR, output_bai)
            output_bam_filename = os.path.join(output_dir.path, output_bam)
            output_bai_filename = os.path.join(output_dir.path, output_bai)
            connor.main(["program_name",
                          input_bam_filename,
                          output_bam_filename])
            self.assertTrue(filecmp.cmp(expect_bam_filename,
                                        output_bam_filename),
                            "{} does not match expected".format(output_bam))
            self.assertTrue(filecmp.cmp(expect_bai_filename,
                                        output_bai_filename),
                            "{} does not match expected".format(output_bai))
