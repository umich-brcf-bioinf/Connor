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
try:
    from builtins import range as iter_range
except ImportError:
    from __builtin__ import xrange as iter_range
from collections import defaultdict, Counter
from copy import deepcopy
from functools import partial
import operator
import sys
import traceback
import time

from sortedcontainers import SortedSet

import connor
import connor.command_validator as command_validator
import connor.samtools as samtools
import connor.consam.pysamwrapper as pysamwrapper
import connor.familyhandler as familyhandler
import connor.utils as utils
from connor.samtools import LoggingWriter

__version__ = connor.__version__

DEFAULT_CONSENSUS_FREQ_THRESHOLD=0.6
DEFAULT_MIN_FAMILY_SIZE_THRESHOLD = 3
DEFAULT_UMT_DISTANCE_THRESHOLD = 1


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


class _TagFamily(object):
    umi_sequence = 0

    def __init__(self,
                 umt,
                 alignments,
                 inexact_match_count,
                 consensus_threshold,
                 family_filter=lambda x: None):
        self.umi_sequence = _TagFamily.umi_sequence
        _TagFamily.umi_sequence += 1
        self._umt = umt
        (self.distinct_cigar_count,
         majority_cigar) = _TagFamily._get_dominant_cigar_stats(alignments)
        self.align_pairs = alignments
        self._mark_minority_cigar(majority_cigar)
        self.inexact_match_count = inexact_match_count
        self.consensus_threshold = consensus_threshold
        self.consensus = self._build_consensus(umt, self.align_pairs)
        self.included_pair_count = sum([1 for p in self.align_pairs if not p.filter_value])
        self.filter_value = family_filter(self)

    def umt(self, format_string=None):
        if format_string:
            return format_string.format(left=self._umt[0], right=self._umt[1])
        else:
            return self._umt

    @staticmethod
    def _get_cigarstring_tuple(paired_alignment):
        return (paired_alignment.left.cigarstring,
                paired_alignment.right.cigarstring)

    def _mark_minority_cigar(self, majority_cigar):
        for pair in self.align_pairs:
            if _TagFamily._get_cigarstring_tuple(pair) != majority_cigar:
                pair.left.filter_value = "minority CIGAR"
                pair.right.filter_value = "minority CIGAR"

    def _generate_consensus_sequence(self, alignment_pairs):
        def consensus_sequence(alignments):
            consensus_seq = self._simple_consensus_sequence(alignments)
            if not consensus_seq:
                consensus_seq = self._complex_consensus_sequence(alignments)
            return consensus_seq
        left_alignments = []
        right_alignments = []
        for align_pair in alignment_pairs:
            left_alignments.append(align_pair.left)
            right_alignments.append(align_pair.right)
        left_consensus_seq = consensus_sequence(left_alignments)
        right_consensus_seq = consensus_sequence(right_alignments)
        return (left_consensus_seq, right_consensus_seq)

    def _simple_consensus_sequence(self, alignments):
        seq_counter = Counter([a.query_sequence for a in alignments])
        consensus_seq = None
        if len(seq_counter) == 1:
            consensus_seq = next(iter(alignments)).query_sequence
        elif len(seq_counter) == 2:
            majority_seq, majority_count = seq_counter.most_common(1)[0]
            total = len(alignments)
            if majority_count / total > self.consensus_threshold:
                consensus_seq = majority_seq
        return consensus_seq

    def _complex_consensus_sequence(self, alignments):
        consensus = []
        for i in iter_range(0, len(alignments[0].query_sequence)):
            counter = Counter([s.query_sequence[i:i+1] for s in alignments])
            base = counter.most_common(1)[0][0]
            freq = counter[base] / sum(counter.values())
            if freq >= self.consensus_threshold:
                consensus.append(base)
            else:
                consensus.append("N")
        return "".join(consensus)

    @staticmethod
    def _select_template_alignment_pair(alignment_pairs):
        top_alignment_pair = None
        best_template = (0, None)
        for alignment_pair in alignment_pairs:
            query_name = alignment_pair.left.query_name
            qual_sum = alignment_pair.left.mapping_quality + \
                    alignment_pair.right.mapping_quality
            if (-qual_sum, query_name) < best_template:
                best_template = (-qual_sum, query_name)
                top_alignment_pair = alignment_pair
        return top_alignment_pair


    @staticmethod
    def _get_dominant_cigar_stats(alignments):
        counter = Counter([_TagFamily._get_cigarstring_tuple(s) for s in alignments])
        number_distict_cigars = len(counter)
        top_two_cigar_count = counter.most_common(2)
        dominant_cigar = top_two_cigar_count[0][0]
        dominant_cigar_count = top_two_cigar_count[0][1]
        if number_distict_cigars > 1 and dominant_cigar_count == top_two_cigar_count[1][1]:
            dominant_cigar = sorted(counter.most_common(),
                                    key=lambda x: (-x[1], x[0]))[0][0]
        return number_distict_cigars, dominant_cigar

    def is_consensus_template(self, connor_align):
        return self.consensus.left.query_name == connor_align.query_name

    def _build_consensus(self, umt, align_pairs):
        included_pairs = [p for p in align_pairs if not p.filter_value]
        template_pair = _TagFamily._select_template_alignment_pair(included_pairs)

        left_align = deepcopy(template_pair.left, {})
        right_align = deepcopy(template_pair.right, {})
        (left_sequence,
         right_sequence) = self._generate_consensus_sequence(included_pairs)
        left_align.query_sequence = left_sequence
        right_align.query_sequence = right_sequence
        left_align.query_qualities = \
                template_pair.left.query_qualities
        right_align.query_qualities = \
                template_pair.right.query_qualities
        consensus_pair = samtools.PairedAlignment(left_align,
                                                   right_align,
                                                   tag_length=len(umt[0]))
        consensus_pair.replace_umt(umt)
        return consensus_pair

#TODO: cgates: reduce complexity
def _build_coordinate_pairs(connor_alignments, excluded_writer):
    MISSING_MATE_FILTER = 'read mate was missing or excluded'
    coords = defaultdict(dict)
    for alignment in connor_alignments:
        if alignment.orientation == 'left':
            key = (alignment.reference_id, alignment.next_reference_start)
            coords[key][alignment.query_name] = alignment
        elif alignment.orientation == 'neither':
            key = (alignment.reference_id, alignment.next_reference_start)
            if key in coords and alignment.query_name in coords[key]:
                align1 = coords[key].pop(alignment.query_name)
                yield samtools.PairedAlignment(align1, alignment)
            else:
                coords[key][alignment.query_name] = alignment
        else:
            key = (alignment.reference_id, alignment.reference_start)
            coord = coords[key]
            l_align = coord.pop(alignment.query_name, None)
            # Clear empty coordinate dict
            if not len(coord):
                del coords[key]
            if l_align:
                yield samtools.PairedAlignment(l_align, alignment)
            else:
                alignment.filter_value = MISSING_MATE_FILTER
                excluded_writer.write(None, None, alignment)
    for aligns in coords.values():
        for align in aligns.values():
            align.filter_value = MISSING_MATE_FILTER
            excluded_writer.write(None, None, align)

class _CoordinateFamilyHolder(object):
    '''Encapsulates how stream of paired aligns are iteratively released as
    sets of pairs which share the same coordinate (coordinate families)'''
    def __init__(self):
        self._coordinate_family = defaultdict(partial(defaultdict, list))
        self._right_coords_in_progress = defaultdict(SortedSet)
        self.pending_pair_count = 0
        self.pending_pair_peak_count = 0

    def _add(self, pair):
        def _start(align):
            return (align.reference_name, align.reference_start)
        self._right_coords_in_progress[pair.right.reference_name].add(pair.right.reference_start)
        right_coord = self._coordinate_family[_start(pair.right)]
        right_coord[_start(pair.left)].append(pair)
        self.pending_pair_count += 1
        self.pending_pair_peak_count = max(self.pending_pair_count,
                                           self.pending_pair_peak_count)

    def _completed_families(self, reference_name, rightmost_boundary):
        '''returns one or more families whose end < rightmost boundary'''
        in_progress = self._right_coords_in_progress[reference_name]
        while len(in_progress):
            right_coord = in_progress[0]
            if right_coord < rightmost_boundary:
                in_progress.pop(0)
                left_families = self._coordinate_family.pop((reference_name, right_coord), {})
                for family in sorted(left_families.values(),
                                     key=lambda x:x[0].left.reference_start):
                    family.sort(key=lambda x: x.query_name)
                    self.pending_pair_count -= len(family)
                    yield family
            else:
                break

    def _remaining_families(self):
        for left_families in self._coordinate_family.values():
            for family in left_families.values():
                self.pending_pair_count -= len(family)
                yield family
            left_families.clear()
        self._coordinate_family.clear()

    #TODO: cgates: reduce complexity
    def build_coordinate_families(self, paired_aligns):
        '''Given a stream of paired aligns, return a list of pairs that share
        same coordinates (coordinate family).  Flushes families in progress
        when any of:
        a) incoming right start > family end
        b) incoming chrom != current chrom
        c) incoming align stream is exhausted'''
        rightmost_start = None
        current_chrom = None
        def _new_coordinate(pair):
            return pair.right.reference_start != rightmost_start
        def _new_chrom(pair):
            return current_chrom != pair.right.reference_name

        for pair in paired_aligns:
            if rightmost_start is None:
                rightmost_start = pair.right.reference_start
                current_chrom = pair.right.reference_name
            if _new_chrom(pair):
                self._right_coords_in_progress[current_chrom].clear()
                rightmost_start = None
                current_chrom = None
                for family in self._remaining_families():
                    yield family
            elif _new_coordinate(pair):
                right = pair.right
                for family in self._completed_families(right.reference_name,
                                                       right.reference_start):
                    yield family
            self._add(pair)

        for family in self._remaining_families():
            yield family

def _build_tag_families(tagged_paired_aligns,
                        ranked_tags,
                        hamming_threshold,
                        consensus_threshold,
                        family_filter=lambda x: None):
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
        tag_family = _TagFamily(tag,
                               tag_aligns[tag],
                               tag_inexact_match_count[tag],
                               consensus_threshold,
                               family_filter)
        tag_families.append(tag_family)
    return tag_families

def _hamming_dist(str1, str2):
    assert len(str1) == len(str2)
    return sum(utils.iter_map(operator.ne, str1, str2))

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
    parser.add_argument('--simplify_pg_header',
                        action="store_true",
                        default=False,
                        help=argparse.SUPPRESS)
    args = parser.parse_args(arguments[1:])
    args.original_command_line = arguments
    if not args.log_file:
        args.log_file = args.output_bam + ".log"
    return args

def _rank_tags(tagged_paired_aligns):
    '''Return the list of tags ranked from most to least popular.'''
    tag_count_dict = defaultdict(int)
    for paired_align in tagged_paired_aligns:
        tag_count_dict[paired_align.umt] += 1
    tags_by_count = utils.sort_dict(tag_count_dict)
    ranked_tags = [tag_count[0] for tag_count in tags_by_count]
    return ranked_tags

def _progress_logger(base_generator,
                     total_rows,
                     log,
                     supplemental_log=lambda x: None):
    row_count = 0
    next_breakpoint = 0
    for item in base_generator:
        row_count += 1
        progress = 100 * row_count / total_rows
        if progress >= next_breakpoint and progress < 100:
            log.info("{}% ({}/{}) alignments processed",
                     next_breakpoint,
                     row_count,
                     total_rows)
            supplemental_log(log)
            next_breakpoint = 10 * int(progress/10) + 10
        yield item
    log.info("100% ({}/{}) alignments processed", row_count, total_rows)
    supplemental_log(log)

def _build_family_filter(args):
    min_family_size = args.min_family_size_threshold
    too_small_msg = 'family too small (<{})'.format(min_family_size)
    def family_size_filter(family):
        if family.included_pair_count < min_family_size:
            return too_small_msg
        else:
            return None
    return family_size_filter

def _build_supplemental_log(coordinate_holder):
    def supplemental_progress_log(log):
        log.debug("{}mb peak memory", utils.peak_memory())
        log.debug("{} pending alignment pairs; {} peak pairs",
                  coordinate_holder.pending_pair_count,
                  coordinate_holder.pending_pair_peak_count)
    return supplemental_progress_log

def _dedup_alignments(args, consensus_writer, annotated_writer, log):
    log.info('reading input bam [{}]', args.input_bam)
    total_aligns = samtools.total_align_count(args.input_bam)
    family_filter = _build_family_filter(args)
    handlers = familyhandler.build_family_handlers(args,
                                                   consensus_writer,
                                                   annotated_writer,
                                                   log)

    bamfile = pysamwrapper.alignment_file(args.input_bam, 'rb')
    coord_family_holder = _CoordinateFamilyHolder()
    supplemental_log = _build_supplemental_log(coord_family_holder)
    progress_gen = _progress_logger(bamfile.fetch(),
                                    total_aligns,
                                    log,
                                    supplemental_log)
    filtered_aligns_gen = samtools.filter_alignments(progress_gen,
                                                     annotated_writer)
    paired_align_gen = _build_coordinate_pairs(filtered_aligns_gen,
                                                    annotated_writer)
    coord_family_gen = coord_family_holder.build_coordinate_families(paired_align_gen)
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
        bam_tags = samtools._build_bam_tags()
        base_annotated_writer = samtools.build_writer(args.input_bam,
                                                      args.annotated_output_bam,
                                                      bam_tags,
                                                      args)
        annotated_writer = LoggingWriter(base_annotated_writer, log)
        consensus_writer = samtools.build_writer(args.input_bam,
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
