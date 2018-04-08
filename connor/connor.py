#! /usr/bin/env python
'''Deduplicates BAM file based on custom inline DNA barcodes.
Emits a new BAM file reduced to a single consensus read for each family of
original reads.'''
##   Copyright 2014 Bioinformatics Core, University of Michigan
##
##   Licensed under the Apache License, Version 2.0 (the "License");
##   you may not use this file except in compliance with the License.
##   You may obtain a copy of the License at
##
##       http://www.apache.org/licenses/LICENSE-2.0
##
##   Unless required by applicable law or agreed to in writing, software
##   distributed under the License is distributed on an "AS IS" BASIS,
##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##   See the License for the specific language governing permissions and
##   limitations under the License.
from __future__ import print_function, absolute_import, division
import argparse
from collections import defaultdict
from functools import partial
import operator
import sys
import traceback
import time

import connor
import connor.command_validator as command_validator
import connor.consam.bamtag as bamtag
import connor.consam.pysamwrapper as pysamwrapper
import connor.consam.readers as readers
import connor.consam.writers as writers
from connor.family import CoordinateFamilyHolder
from connor.family import TagFamily
import connor.familyhandler as familyhandler
import connor.utils as utils

__version__ = connor.__version__

DEFAULT_CONSENSUS_FREQ_THRESHOLD=0.6
DEFAULT_MIN_FAMILY_SIZE_THRESHOLD = 3
DEFAULT_UMT_DISTANCE_THRESHOLD = 1
FILTER_ORDERS = ['count', 'name']
DEFAULT_FILTER_ORDER = FILTER_ORDERS[0]

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
        raise utils.UsageError(message)


def _build_tag_families(tagged_paired_aligns,
                        ranked_tags,
                        hamming_threshold,
                        consensus_threshold,
                        family_filter=lambda _: None):
    '''Partition paired aligns into families.

    Each read is considered against each ranked tag until all reads are
    partitioned into families.'''
    tag_aligns = defaultdict(set)
    tag_inexact_match_count = defaultdict(int)

    for paired_align in tagged_paired_aligns:
        (left_umt, right_umt) =  paired_align.umt
        for best_tag in ranked_tags:
            if paired_align.umt == best_tag:
                tag_aligns[best_tag].add(paired_align)
                break
            elif left_umt == best_tag[0] or right_umt == best_tag[1]:
                tag_aligns[best_tag].add(paired_align)
                tag_inexact_match_count[best_tag] += 1
                break
            elif (_hamming_dist(left_umt, best_tag[0]) <= hamming_threshold) \
                or (_hamming_dist(right_umt, best_tag[1]) <= hamming_threshold):
                tag_aligns[best_tag].add(paired_align)
                tag_inexact_match_count[best_tag] += 1
                break
    tag_families = []
    for tag in sorted(tag_aligns):
        tag_family = TagFamily(tag,
                               tag_aligns[tag],
                               tag_inexact_match_count[tag],
                               consensus_threshold,
                               family_filter)
        tag_families.append(tag_family)
    return tag_families

def _hamming_dist(str1, str2):
    assert len(str1) == len(str2)
    return sum(utils.iter_map(operator.ne, str1, str2))

def _rank_tags(tagged_paired_aligns):
    '''Return the list of tags ranked from most to least popular.'''
    tag_count_dict = defaultdict(int)
    for paired_align in tagged_paired_aligns:
        tag_count_dict[paired_align.umt] += 1
    tags_by_count = utils.sort_dict(tag_count_dict)
    ranked_tags = [tag_count[0] for tag_count in tags_by_count]
    return ranked_tags

def _build_family_filter(args):
    min_family_size = args.min_family_size_threshold
    too_small_msg = 'family too small (<{})'.format(min_family_size)
    def _family_size_filter(family):
        if family.included_pair_count < min_family_size:
            return too_small_msg
        return None
    return _family_size_filter

def _build_supplemental_log(coordinate_holder):
    def _supplemental_progress_log(log):
        log.debug("{}mb peak memory", utils.peak_memory())
        log.debug("{} pending alignment pairs; {} peak pairs",
                  coordinate_holder.pending_pair_count,
                  coordinate_holder.pending_pair_peak_count)
    return _supplemental_progress_log

#TODO: cgates push bamfile handling into readers.paired_reder
def _dedup_alignments(args, consensus_writer, annotated_writer, log):
    log.info('reading input bam [{}]', args.input_bam)

    bamfile = pysamwrapper.alignment_file(args.input_bam, 'rb')
    total_aligns = pysamwrapper.total_align_count(args.input_bam)
    coord_family_holder = CoordinateFamilyHolder()
    supplemental_log = _build_supplemental_log(coord_family_holder)
    paired_align_gen = readers.paired_reader(bamfile,
                                             total_aligns,
                                             log,
                                             supplemental_log,
                                             annotated_writer)
    coord_family_gen = coord_family_holder.build_coordinate_families(paired_align_gen)

    family_filter = _build_family_filter(args)
    handlers = familyhandler.build_family_handlers(args,
                                                   consensus_writer,
                                                   annotated_writer,
                                                   log)
    for coord_family in coord_family_gen:
        ranked_tags = _rank_tags(coord_family)
        tag_families = _build_tag_families(coord_family,
                                           ranked_tags,
                                           args.umt_distance_threshold,
                                           args.consensus_freq_threshold,
                                           family_filter)
        for handler in handlers:
            for tag_family in tag_families:
                handler.handle(tag_family)

    for handler in handlers:
        handler.end()

    bamfile.close()

def _parse_command_line_args(arguments):
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
    parser.add_argument('--simplify_pg_header',
                        action="store_true",
                        default=False,
                        help=argparse.SUPPRESS)
    args = parser.parse_args(arguments[1:])
    args.original_command_line = arguments
    if not args.log_file:
        args.log_file = args.output_bam + ".log"
    return args


def main(command_line_args=None):
    '''Connor entry point.  See help for more info'''
    log = None
    if not command_line_args:
        command_line_args = sys.argv
    try:
        start_time = time.time()
        args = _parse_command_line_args(command_line_args)
        log = utils.Logger(args)
        command_validator.preflight(args, log)
        log.info('connor begins (v{})', __version__)
        log.info('logging to [{}]', args.log_file)
        utils.log_environment_info(log, args)
        bam_tags = bamtag.build_bam_tags()
        base_annotated_writer = writers.build_writer(args.input_bam,
                                                      args.annotated_output_bam,
                                                      bam_tags,
                                                      args)
        sort_filter_by_name = args.filter_order == 'name'
        annotated_writer = writers.LoggingWriter(base_annotated_writer,
                                                 log,
                                                 sort_filter_by_name)
        consensus_writer = writers.build_writer(args.input_bam,
                                                 args.output_bam,
                                                 bam_tags,
                                                 args)
        _dedup_alignments(args, consensus_writer, annotated_writer, log)
        annotated_writer.close(log)
        consensus_writer.close(log)
        warning = ' **See warnings above**' if log.warning_occurred else ''
        elapsed_time = int(time.time() - start_time)
        log.info("connor complete ({} seconds, {}mb peak memory).{}",
             elapsed_time,
             utils.peak_memory(),
             warning)
    except utils.UsageError as usage_error:
        message = "connor usage problem: {}".format(str(usage_error))
        print(message, file=sys.stderr)
        print("See 'connor --help'.", file=sys.stderr)
        sys.exit(1)
    except Exception: #pylint: disable=broad-except
        if log:
            show = log.error
        else:
            show = partial(print, file=sys.stderr)
        show("An unexpected error occurred")
        show(traceback.format_exc())
        exit(1)


if __name__ == '__main__':
    main(sys.argv)
