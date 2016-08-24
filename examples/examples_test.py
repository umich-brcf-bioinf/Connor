from __future__ import print_function, absolute_import
import filecmp
import os
from test.utils_test import BaseConnorTestCase

from testfixtures import TempDirectory

import connor.connor as connor

INPUT_DIR=os.path.realpath(os.path.dirname(__file__))

class ExamplesFunctionalTest(BaseConnorTestCase):
    def test_examples(self):
        with TempDirectory() as output_dir:
            def out_path(filename):
                return os.path.join(output_dir.path, filename)
            def in_path(filename):
                return os.path.join(INPUT_DIR, filename)
            input_bam = "PIK3CA-original.bam"
            output_deduped_bam = "PIK3CA-deduped.bam"
            output_deduped_bai = output_deduped_bam + ".bai"
            output_annotated_bam = "PIK3CA-annotated.bam"
            output_annotated_bai = output_annotated_bam + ".bai"

            input_bam_path = in_path(input_bam)
            expect_deduped_bam_path = in_path(output_deduped_bam)
            expect_deduped_bai_path = in_path(output_deduped_bai)
            expect_annotated_bam_path = in_path(output_annotated_bam)
            expect_annotated_bai_path = in_path(output_annotated_bai)

            output_deduped_bam_path = out_path(output_deduped_bam)
            output_deduped_bai_path = out_path(output_deduped_bai)
            output_annotated_bam_path = out_path(output_annotated_bam)
            output_annotated_bai_path = out_path(output_annotated_bai)

            connor.main(["program_name",
                          input_bam_path,
                          output_deduped_bam_path,
                          "--annotated_output_bam",
                          output_annotated_bam_path,
                          "--verbose",
                          "--simplify_pg_header"])

            def compare_files(expected_path, output_path):
                expected_name = os.path.basename(expected_path)
                output_name = os.path.basename(output_path)
                if not os.path.isfile(output_path):
                    result = 'output file {} not found'.format(output_name)
                elif not os.path.isfile(expected_path):
                    result = 'expected file {} not found'.format(expected_name)
                elif not filecmp.cmp(expected_path, output_path, shallow=False):
                    msg = '{} does not match expected'
                    result = msg.format(output_path)
                else:
                    result = 'OK'
                return expected_name, result

            comparisons = [(expect_deduped_bam_path, output_deduped_bam_path),
                           (expect_deduped_bai_path, output_deduped_bai_path),
                           (expect_annotated_bam_path, output_annotated_bam_path),
                           (expect_annotated_bai_path, output_annotated_bai_path)]
            diff_files = dict([compare_files(expected, output) for expected, output in comparisons])

        files_match = set(diff_files.values()) == set(["OK"])
        self.assertTrue(files_match, diff_files)
