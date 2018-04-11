'''Parses command line args'''
from __future__ import print_function, absolute_import, division
import argparse

from connor import __version__
from connor.utils import UsageError

FILTER_ORDERS = ['count', 'name']

DEFAULT_CONSENSUS_FREQ_THRESHOLD=0.6
DEFAULT_FILTER_ORDER = FILTER_ORDERS[0]
DEFAULT_MIN_FAMILY_SIZE_THRESHOLD = 3
DEFAULT_UMT_DISTANCE_THRESHOLD = 1
DEFAULT_UMT_LENGTH = 6

DESCRIPTION=\
'''Deduplicates BAM file based on custom inline DNA barcodes.
Emits a new BAM file reduced to a single consensus read for each family of
original reads.
'''

class _ConnorArgumentParser(argparse.ArgumentParser):
    '''Argument parser that raises UsageError instead of exiting.'''
    #pylint: disable=too-few-public-methods
    def error(self, message):
        '''Suppress default exit behavior'''
        raise UsageError(message)


def parse_command_line_args(arguments):
    parser = _ConnorArgumentParser( \
        formatter_class=argparse.RawTextHelpFormatter,
        usage="connor input_bam output_bam",
        description=(DESCRIPTION),
        epilog="v"+__version__)

    parser.add_argument("-V",
                        "--version",
                        action='version',
                        version=__version__)
    parser.add_argument("-v",
                        "--verbose",
                        action="store_true",
                        default=False,
                        help="print all log messages to console")
    parser.add_argument('input_bam',
                        help="path to input BAM")
    parser.add_argument('output_bam',
                        help="path to deduplicated output BAM")
    parser.add_argument('--force',
                        action="store_true",
                        default=False,
                        help="=False. Override validation warnings")
    parser.add_argument('--log_file',
                        type=str,
                        help="={output_filename}.log. Path to verbose log file")
    parser.add_argument('--annotated_output_bam',
                        type=str,
                        help=("path to output BAM containing all original "
                              "aligns annotated with BAM tags"))
    parser.add_argument("-f", "--consensus_freq_threshold",
                        type=float,
                        default = DEFAULT_CONSENSUS_FREQ_THRESHOLD,
                        help = \
"""={} (0..1.0): Ambiguous base calls at a specific position in a family are
 transformed to either majority base call, or N if the majority percentage
 is below this threshold. (Higher threshold results in more Ns in
 consensus.)""".format(DEFAULT_CONSENSUS_FREQ_THRESHOLD))
    parser.add_argument("-s", "--min_family_size_threshold",
                        type=int,
                        default = DEFAULT_MIN_FAMILY_SIZE_THRESHOLD,
                        help=\
"""={} (>=0): families with count of original reads < threshold are excluded
 from the deduplicated output. (Higher threshold is more
 stringent.)""".format(DEFAULT_MIN_FAMILY_SIZE_THRESHOLD))
    parser.add_argument("-d", "--umt_distance_threshold",
                        type=int,
                        default = DEFAULT_UMT_DISTANCE_THRESHOLD,
                        help=\
"""={} (>=0); UMTs equal to or closer than this Hamming distance will be
 combined into a single family. Lower threshold make more families with more
 consistent UMTs; 0 implies UMI must match
 exactly.""".format(DEFAULT_UMT_DISTANCE_THRESHOLD))
    parser.add_argument("--filter_order",
                        choices=FILTER_ORDERS,
                        default = DEFAULT_FILTER_ORDER,
                        help=\
"""={}; determines how filters will be ordered in the log
 results""".format(DEFAULT_FILTER_ORDER))
    parser.add_argument("--umt_length",
                        type=int,
                        default = DEFAULT_UMT_LENGTH,
                        help=\
"""={} (>=1); length of UMT""".format(DEFAULT_UMT_LENGTH))
    parser.add_argument('--simplify_pg_header',
                        action="store_true",
                        default=False,
                        help=argparse.SUPPRESS)
    args = parser.parse_args(arguments[1:])
    args.original_command_line = arguments
    if not args.log_file:
        args.log_file = args.output_bam + ".log"
    return args
