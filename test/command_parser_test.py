#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments
#pylint: disable=deprecated-method
from __future__ import print_function, absolute_import, division
from test.utils_test import BaseConnorTestCase

import connor.command_parser as command_parser
from connor.utils import UsageError
from connor.command_parser import parse_command_line_args

class CommandParser(BaseConnorTestCase):
    def test_parse_command_line_args(self):
        namespace = parse_command_line_args(["command",
                                             "input.bam",
                                             "output.bam"])
        self.assertEqual("input.bam", namespace.input_bam)
        self.assertEqual("output.bam", namespace.output_bam)
        self.assertEqual(False, namespace.force)
        self.assertEqual(False, namespace.simplify_pg_header)
        self.assertEqual(False, namespace.verbose)
        self.assertEqual("output.bam.log", namespace.log_file)
        self.assertEqual(None, namespace.annotated_output_bam)
        self.assertEqual(command_parser.DEFAULT_CONSENSUS_FREQ_THRESHOLD,
                         namespace.consensus_freq_threshold)
        self.assertEqual(command_parser.DEFAULT_MIN_FAMILY_SIZE_THRESHOLD,
                         namespace.min_family_size_threshold)
        self.assertEqual(command_parser.DEFAULT_UMT_DISTANCE_THRESHOLD,
                         namespace.umt_distance_threshold)
        self.assertEquals(['command', 'input.bam', 'output.bam'],
                          namespace.original_command_line)
        self.assertEquals(command_parser.DEFAULT_FILTER_ORDER,
                          namespace.filter_order)
        self.assertEquals(command_parser.DEFAULT_UMT_LENGTH,
                          namespace.umt_length)
        self.assertEqual(13, len(vars(namespace)))

    def test_parse_command_line_args_throwsConnorUsageError(self):
        self.assertRaises(UsageError,
                          parse_command_line_args,
                          ["command", "input"])
        self.assertRaises(UsageError,
                          parse_command_line_args,
                          ["command",
                           "input",
                           "output",
                           "something else"])
