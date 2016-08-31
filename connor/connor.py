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
import os
import platform
import sys
import traceback
import resource
import time

import pysam
from sortedcontainers import SortedSet

import connor
import connor.command_validator as command_validator
import connor.samtools as samtools
import connor.familyhandler as familyhandler
import connor.utils as utils
from connor.samtools import LoggingWriter

__version__ = connor.__version__

DEFAULT_TAG_LENGTH = 6
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


class _LightweightAlignment(object):
    '''Minimal info from PySam.AlignedSegment used to expedite pos grouping.'''
    def __init__(self, aligned_segment):
        self.name = aligned_segment.query_name
        chrom = aligned_segment.reference_name
        pos1 = aligned_segment.reference_start
        pos2 = aligned_segment.next_reference_start
        self.reference_end = aligned_segment.reference_end
        if pos1 < pos2:
            self.key = (chrom, pos1, pos2)
            self.left_pos = pos1
        else:
            self.key = (chrom, pos2, pos1)
            self.left_pos = pos2


class _LightweightPair(object):
    '''Minimal info from PySam.AlignedSegment used to expedite pos grouping.'''
    def __init__(self, aligned_segment1, aligned_segment2):
        self.name = aligned_segment1.query_name
        chrom = aligned_segment1.reference_name
        left_start = min(aligned_segment1.reference_start,
                         aligned_segment2.reference_start)
        right_end = max(aligned_segment1.reference_end,
                        aligned_segment2.reference_end)
        self.key = (chrom, left_start, right_end)


class _PairedAlignment(object):
    '''Represents the left and right align pairs from an single sequence.'''
    def __init__(self,
                 left_alignment,
                 right_alignment,
                 tag_length=DEFAULT_TAG_LENGTH):
        if left_alignment.query_name != right_alignment.query_name:
            msg = 'Inconsistent query names ({} != {})'
            raise ValueError(msg.format(left_alignment.query_name,
                                        right_alignment.query_name))
        self.query_name = left_alignment.query_name
        self.left = left_alignment
        self.right = right_alignment
        self._tag_length = tag_length
        left_umt = self.left.query_sequence[0:self._tag_length]
        right_umt = self.right.query_sequence[-1 * self._tag_length:]
        self.umt = (left_umt, right_umt)

    @property
    def filter_value(self):
        if self.left.filter_value or self.right.filter_value:
            return (self.left.filter_value, self.right.filter_value)
        else:
            return None

    def replace_umt(self, umt):
        def _byte_array_to_string(sequence):
            if isinstance(sequence, str):
                return sequence
            else:
                return str(sequence.decode("utf-8"))
        if not (umt[0] or umt[1]) or \
            (len(umt[0]) != self._tag_length) or \
            (len(umt[1]) != self._tag_length):
            msg = "Each UMT must match tag_length ({})"
            raise ValueError(msg.format(self._tag_length))
        left_qual = self.left.query_qualities
        right_qual = self.right.query_qualities
        left_query_frag = self.left.query_sequence[len(umt[0]):]
        left_query_frag_str = _byte_array_to_string(left_query_frag)
        self.left.query_sequence = umt[0] + left_query_frag_str
        right_query_frag = self.right.query_sequence[:-len(umt[1])]
        right_query_frag_str = _byte_array_to_string(right_query_frag)
        self.right.query_sequence = right_query_frag_str + umt[1]
        self.umt = umt
        self.left.query_qualities = left_qual
        self.right.query_qualities = right_qual

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):
        return hash(self.left) * hash(self.right)

    def __repr__(self):
        return ("Pair({}|{}|{}, "
                "{}|{}|{})").format(self.left.query_name,
                                    self.left.reference_start,
                                    self.left.query_sequence,
                                    self.right.query_name,
                                    self.right.reference_start,
                                    self.right.query_sequence)

#TODO: cgates: move consensus builder into separate class/module
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
        self.umt = umt
        (self.distinct_cigar_count,
         majority_cigar) = _TagFamily._get_dominant_cigar_stats(alignments)
        self.align_pairs = alignments
        self._mark_minority_cigar(majority_cigar)
        self.inexact_match_count = inexact_match_count
        self.consensus_threshold = consensus_threshold
        self.consensus = self._build_consensus(umt, self.align_pairs)
        self.included_pair_count = sum([1 for p in self.align_pairs if not p.filter_value])
        self.filter_value = family_filter(self)

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


    #TODO: (cgates) tags should not assume umt is a tuple and symmetric
    #between left and right
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
        consensus_pair = _PairedAlignment(left_align,
                                          right_align,
                                          tag_length=len(umt[0]))
        consensus_pair.replace_umt(umt)
        return consensus_pair


def _build_coordinate_read_name_manifest(lw_aligns):
    '''Return a dict mapping coordinates to set of aligned querynames.

    Constructed on a preliminary pass through the input BAM, this lightweight
    dict informs downstream processing that the collection of reads at a
    coordinate can be released.
    '''
    af_dict = defaultdict(set)
    for lwa in lw_aligns:
        af_dict[lwa.key].add(lwa.name)
    return af_dict

def _build_coordinate_families(aligned_segments,
                               coord_read_name_manifest,
                               excluded_writer):
    '''Generate sets of PairedAlignments that share the same coordinates.'''
    family_dict = defaultdict(set)
    pairing_dict = {}
    for aseg in aligned_segments:
        if not aseg.query_name in pairing_dict:
            pairing_dict[aseg.query_name]= aseg
        else:
            paired_align = _PairedAlignment(pairing_dict.pop(aseg.query_name),
                                            aseg)
            key = _LightweightPair(paired_align.left,
                                   paired_align.right).key
            family_dict[key].add(paired_align)
            coord_read_name_manifest[key].remove(aseg.query_name)
            if not coord_read_name_manifest[key]:
                yield family_dict.pop(key)

    for align in sorted(pairing_dict.values(), key=lambda a:a.query_name):
        align.filter_value = 'read mate was missing or excluded'
        excluded_writer.write(None, align)

def _build_coordinate_pairs_deux(connor_alignments):
    coords = defaultdict(dict)
    for alignment in connor_alignments:
        if alignment.orientation == 'left':
            key = (alignment.reference_id, alignment.next_reference_start)
            coords[key][alignment.query_name] = alignment
        else:
            key = (alignment.reference_id, alignment.reference_start)
            #TODO: cgates: this line will raise on orphaned right read
            l_align = coords[key].pop(alignment.query_name)
            yield _PairedAlignment(l_align, alignment)
    #TODO: cgates: write orphaned pairs to excluded writer

class _CoordinateFamilyHolder(object):
    '''Encapsulates how stream of paired aligns are iteratively released as
    sets of pairs which share the same coordinate (coordinate families)'''
    def __init__(self):
        self.coordinate_family = defaultdict(partial(defaultdict, list))
        self.right_coords_in_progress = SortedSet()

    def _add(self, pair):
        right = pair.right.reference_end
        left = pair.left.reference_start
        self.right_coords_in_progress.add(right)
        self.coordinate_family[right][left].append(pair)

    def _completed_families(self, rightmost_boundary):
        '''returns one or more families whose end < rightmost boundary'''
        while len(self.right_coords_in_progress):
            right_coord = self.right_coords_in_progress[0]
            if right_coord < rightmost_boundary:
                self.right_coords_in_progress.pop(0)
                left_families = self.coordinate_family.pop(right_coord)
                for family in left_families.values():
                    yield family
            else:
                break

    def _remaining_families(self):
        for left_families in self.coordinate_family.values():
            for family in left_families.values():
                yield family
            left_families.clear()
        self.coordinate_family.clear()

    #TODO: cgates: can we reduce the complexity here?
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
                for family in self._remaining_families():
                    yield family
            elif _new_coordinate(pair):
                rightmost_start = pair.right.reference_start
                for family in self._completed_families(rightmost_start):
                    yield family
            self._add(pair)

        for family in self._remaining_families():
            yield family


def _build_coordinate_families_deux(paired_aligns):
    family_holder = _CoordinateFamilyHolder()
    for family in family_holder.build_coordinate_families(paired_aligns):
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


def _build_lightweight_pairs(aligned_segments):
    name_pairs = dict()
    lightweight_pairs = list()
    total_align_count = 0
    for align_segment in aligned_segments:
        total_align_count += 1
        query_name = align_segment.query_name
        if not query_name in name_pairs:
            name_pairs[align_segment.query_name] = align_segment
        else:
            new_pair = _LightweightPair(name_pairs.pop(query_name),
                                        align_segment)
            lightweight_pairs.append(new_pair)
    return lightweight_pairs

#TODO: cgates: improve how this is tested
def _progress_logger(base_generator, total_rows, log):
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
            log.debug("{}mb peak memory", _peak_memory())
            next_breakpoint = 10 * int(progress/10) + 10
        yield item
    log.info("100% ({}/{}) alignments processed", row_count, total_rows)
    log.debug("{}mb peak memory", _peak_memory())

def _build_family_filter(args):
    min_family_size = args.min_family_size_threshold
    too_small_msg = 'family too small (<{})'.format(min_family_size)
    def family_size_filter(family):
        if family.included_pair_count < min_family_size:
            return too_small_msg
        else:
            return None
    return family_size_filter

def _build_manifest(input_bam):
    bamfile = samtools.alignment_file(input_bam, 'rb')
    counter = utils.CountingGenerator()
    included_aligns = samtools.filter_alignments(counter.count(bamfile.fetch()))
    lightweight_pairs = _build_lightweight_pairs(included_aligns)
    bamfile.close()
    coord_manifest = _build_coordinate_read_name_manifest(lightweight_pairs)
    return counter.item_count, coord_manifest

def _dedup_alignments(args, consensus_writer, annotated_writer, log):
    log.info('reading input bam [{}]', args.input_bam)
    (total_aligns,
     coord_manifest) = _build_manifest(args.input_bam)
    family_filter = _build_family_filter(args)
    handlers = familyhandler.build_family_handlers(args,
                                                   consensus_writer,
                                                   annotated_writer,
                                                   log)

    bamfile = samtools.alignment_file(args.input_bam, 'rb')
    progress_gen = _progress_logger(bamfile.fetch(),
                                    total_aligns,
                                    log)
    filtered_aligns_gen = samtools.filter_alignments(progress_gen,
                                                     annotated_writer)
    for coord_family in _build_coordinate_families(filtered_aligns_gen,
                                                   coord_manifest,
                                                   annotated_writer):
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


def _log_environment_info(log, args):
    log.debug('original_command_line|{}',' '.join(args.original_command_line))
    log.debug('command_options|{}', vars(args))
    log.debug('command_cwd|{}', os.getcwd ())
    log.debug('platform_uname|{}', platform.uname())
    log.debug('platform_python_version|{}', platform.python_version())
    log.debug('pysam_version|{}', pysam.__version__)

#TODO: cgates: None checking/cyclomatic complexity could be simplified with UNPLACED/NULL family object
def _build_bam_tags():
    def combine_filters(family, align):
        filter_values = [x.filter_value for x in [family, align] if x and x.filter_value]
        if filter_values:
            return ";".join(filter_values).replace('; ', ';')
        else:
            return None
    tags = [
        samtools.BamTag("X0", "Z",
                        ("filter (why the alignment was excluded)"),
                        combine_filters),
        samtools.BamTag("X1", "i",
                        "unique identifier for this alignment family",
                        lambda fam, align: fam.umi_sequence if fam else None),
        samtools.BamTag("X2", "Z",
                        ("L~R UMT barcodes for this alignment family; because "
                         "of fuzzy matching the family UMT may be distinct "
                         "from the UMT of the original alignment"),
                        lambda fam, align: "{0}~{1}".format(fam.umt[0],
                                                            fam.umt[1]) if fam else None),
        samtools.BamTag("X3", "i",
                        "family size (number of align pairs in this family)",
                        lambda fam, align: fam.included_pair_count if fam else None),
        samtools.BamTag("X4", "i",
                        ("presence of this tag signals that this alignment "
                         "would be the template for the consensus alignment"),
                        lambda fam, align: 1 if fam and fam.consensus.left.query_name == align.query_name else None)]
    return tags


def _peak_memory():
    peak_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_memory_mb = peak_memory/1024
    if sys.platform == 'darwin':
        peak_memory_mb /= 1024
    return int(peak_memory_mb)


def main(command_line_args=None):
    '''Connor entry point.  See help for more info'''
    if not command_line_args:
        command_line_args = sys.argv
    try:
        start_time = time.time()
        args = _parse_command_line_args(command_line_args)
        log = utils.Logger(args)
        command_validator.preflight(args, log)
        log.info('connor begins (v{})', __version__)
        log.info('logging to [{}]', args.log_file)
        _log_environment_info(log, args)
        bam_tags = _build_bam_tags()
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
             _peak_memory(),
             warning)
    except utils.UsageError as usage_error:
        message = "connor usage problem: {}".format(str(usage_error))
        print(message, file=sys.stderr)
        print("See 'connor --help'.", file=sys.stderr)
        sys.exit(1)
    except Exception: #pylint: disable=broad-except
        log.error("An unexpected error occurred")
        log.error(traceback.format_exc())
        exit(1)


if __name__ == '__main__':
    main(sys.argv)
