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
from collections import defaultdict, Counter
from copy import deepcopy
import operator
import os
import platform
import sys
import traceback
import resource
import time

import pysam

import connor
import connor.samtools as samtools
import connor.familyhandler as familyhandler
import connor.utils as utils

__version__ = connor.__version__

DEFAULT_TAG_LENGTH = 6
DEFAULT_CONSENSUS_FREQ_THRESHOLD=0.6
DEFAULT_MIN_FAMILY_SIZE_THRESHOLD = 3
DEFAULT_UMI_DISTANCE_THRESHOLD = 1


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
    def __init__(self, left_alignment,
                 right_alignment,
                 tag_length=DEFAULT_TAG_LENGTH):
        self.left_alignment = left_alignment
        self.right_alignment = right_alignment
        self._tag_length = tag_length
        left_umt = self.left_alignment.query_sequence[0:self._tag_length]
        right_umt = self.right_alignment.query_sequence[-1 * self._tag_length:]
        self.umi = (left_umt, right_umt)

    def replace_umt(self, umt):
        def _byte_array_to_string(sequence):
            if isinstance(sequence, str):
                return sequence
            else:
                return str(sequence.decode("utf-8"))
        if not (umt[0] or umt[1]) or \
            (len(umt[0]) != self._tag_length) or \
            (len(umt[1]) != self._tag_length):
            raise ValueError("Each UMT must match tag_length ({})".format(self._tag_length))
        left_qual = self.left_alignment.query_qualities
        right_qual = self.right_alignment.query_qualities
        left_query_frag = self.left_alignment.query_sequence[len(umt[0]):]
        left_query_frag_str = _byte_array_to_string(left_query_frag)
        self.left_alignment.query_sequence = umt[0] + left_query_frag_str
        right_query_frag = self.right_alignment.query_sequence[:-len(umt[1])]
        right_query_frag_str = _byte_array_to_string(right_query_frag)
        self.right_alignment.query_sequence = right_query_frag_str + umt[1]
        self.umi = umt
        self.left_alignment.query_qualities = left_qual
        self.right_alignment.query_qualities = right_qual

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):
        return hash(self.left_alignment) * hash(self.right_alignment)

    def __repr__(self):
        return ("Pair({}|{}|{}, "
                "{}|{}|{})").format(self.left_alignment.query_name,
                                    self.left_alignment.reference_start,
                                    self.left_alignment.query_sequence,
                                    self.right_alignment.query_name,
                                    self.right_alignment.reference_start,
                                    self.right_alignment.query_sequence)

#TODO: cgates: move consensus builder into separate class
class _TagFamily(object):
    FILTER_FORMAT = "small family (<{})"
    umi_sequence = 0

    @staticmethod
    def filter_small_alignments(alignments, min_family_size):
        if len(alignments) < min_family_size:
            align_filter = _TagFamily.FILTER_FORMAT.format(min_family_size)
            for align in alignments:
                align.left_alignment.filter = align_filter
                align.right_alignment.filter = align_filter


    def __init__(self,
                 umi,
                 alignments,
                 inexact_match_count,
                 consensus_threshold,
                 min_family_size=DEFAULT_MIN_FAMILY_SIZE_THRESHOLD):
        self.umi_sequence = _TagFamily.umi_sequence
        _TagFamily.umi_sequence += 1
        self.umi = umi
        self.input_alignment_count = len(alignments)
        (self.distinct_cigar_count,
         self.minority_cigar_percentage,
         dominant_cigar) = _TagFamily._get_dominant_cigar_stats(alignments)
        #TODO: cgates: Better to store as single collection of ConnorAligns
        (self.alignments,
         self.excluded_alignments) = _TagFamily._get_alignments_for_dominant_cigar(dominant_cigar,
                                                                                   alignments)
        _TagFamily.filter_small_alignments(self.alignments, min_family_size)
        self.inexact_match_count = inexact_match_count
        #Necessary to make output deterministic
        self.alignments.sort(key=lambda x: x.left_alignment.query_name)
        self.consensus_threshold = consensus_threshold
        self.consensus = self._build_consensus(umi, self.alignments)

    @staticmethod
    def _get_cigarstring_tuple(paired_alignment):
        return (paired_alignment.left_alignment.cigarstring,
                paired_alignment.right_alignment.cigarstring)


    @staticmethod
    def _get_alignments_for_dominant_cigar(dominant_cigar,
                                           list_of_alignments):
        included_alignments = []
        excluded_alignments = []
        for pair in list_of_alignments:
            if _TagFamily._get_cigarstring_tuple(pair) == dominant_cigar:
                included_alignments.append(pair)
            else:
                pair.left_alignment.filter = "minority CIGAR"
                pair.right_alignment.filter = "minority CIGAR"
                excluded_alignments.append(pair)
        return included_alignments, excluded_alignments

#    def _generate_consensus_sequence(self, list_of_alignments):
#        consensus_seq = self._generate_consensus_sequence_quickly(list_of_alignments)
#        if not consensus_seq:
#            consensus_seq = self._generate_consensus_sequence_long(list_of_alignments)
#        return consensus_seq
    
    def _generate_consensus_sequence(self, list_of_alignment_pairs):
        left_alignments = []
        right_alignments = []
        for align_pair in list_of_alignment_pairs:
            left_alignments.append(align_pair.left_alignment)
            right_alignments.append(align_pair.right_alignment)
        left_consensus_seq = self._generate_consensus_sequence_quickly(left_alignments)
        if not left_consensus_seq:
            left_consensus_seq = self._generate_consensus_sequence_long(left_alignments)
        right_consensus_seq = self._generate_consensus_sequence_quickly(right_alignments)
        if not right_consensus_seq:
            right_consensus_seq = self._generate_consensus_sequence_long(right_alignments)
        return (left_consensus_seq, right_consensus_seq)

    def _generate_consensus_sequence_quickly(self, list_of_alignments):
        seq_counter = Counter([a.query_sequence for a in list_of_alignments])
        consensus_seq = None
        if len(seq_counter) == 1:
            consensus_seq = list_of_alignments[0].query_sequence
        elif len(seq_counter) == 2:
            majority_seq, majority_count = seq_counter.most_common(1)[0]
            total = len(list_of_alignments)
            if majority_count / total > self.consensus_threshold:
                consensus_seq = majority_seq
        return consensus_seq

    def _generate_consensus_sequence_long(self, list_of_alignments):
        consensus = []
        for i in utils.zrange(0, len(list_of_alignments[0].query_sequence)):
            counter = Counter([s.query_sequence[i:i+1] for s in list_of_alignments])
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
            query_name = alignment_pair.left_alignment.query_name
            qual_sum = alignment_pair.left_alignment.mapping_quality + \
                    alignment_pair.right_alignment.mapping_quality
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
        if len(top_two_cigar_count) == 1:
            minority_cigar_percentage = 0
        elif dominant_cigar_count == top_two_cigar_count[0][1]:
            dominant_cigar = sorted(counter.most_common(),
                                    key=lambda x: (-x[1], x[0]))[0][0]
            minority_cigar_percentage = top_two_cigar_count[1][1]/len(alignments)
        else:
            minority_cigar_percentage = top_two_cigar_count[1][1]/len(alignments)
        return number_distict_cigars, minority_cigar_percentage, dominant_cigar


    #TODO: (cgates) tags should not assume umi is a tuple and symmetric
    #between left and right
    def _build_consensus(self, umt, align_pairs):
        template_alignment_pair = _TagFamily._select_template_alignment_pair(align_pairs)
        
        left_consensus_align = deepcopy(template_alignment_pair.left_alignment,
                                        {})
        right_consensus_align = deepcopy(template_alignment_pair.right_alignment,
                                         {})
        consensus_align = _PairedAlignment(left_consensus_align,
                                          right_consensus_align,
                                          tag_length=len(umt[0]))
        (left_consensus_sequence,
         right_consensus_sequence) = self._generate_consensus_sequence(align_pairs)
        left_consensus_align.query_sequence = left_consensus_sequence
        right_consensus_align.query_sequence = right_consensus_sequence
        consensus_align.left_alignment.query_qualities = template_alignment_pair.left_alignment.query_qualities
        consensus_align.right_alignment.query_qualities = template_alignment_pair.right_alignment.query_qualities
        consensus_align.replace_umt(umt)
        return consensus_align


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

def _build_coordinate_families(aligned_segments, coord_read_name_manifest):
    '''Generate sets of PairedAlignments that share the same coordinates.'''
    family_dict = defaultdict(set)
    pairing_dict = {}
    for aseg in aligned_segments:
        if not aseg.query_name in pairing_dict:
            pairing_dict[aseg.query_name]= aseg
        else:
            paired_align = _PairedAlignment(pairing_dict.pop(aseg.query_name),
                                            aseg)
            key = _LightweightPair(paired_align.left_alignment,
                                   paired_align.right_alignment).key
            family_dict[key].add(paired_align)
            coord_read_name_manifest[key].remove(aseg.query_name)
            if not coord_read_name_manifest[key]:
                yield family_dict.pop(key)

def _build_tag_families(tagged_paired_aligns,
                        ranked_tags,
                        hamming_threshold,
                        consensus_threshold,
                        min_family_size=DEFAULT_MIN_FAMILY_SIZE_THRESHOLD):
    '''Partition paired aligns into families.

    Each read is considered against each ranked tag until all reads are
    partitioned into families.'''
    tag_aligns = defaultdict(set)
    tag_inexact_match_count = defaultdict(int)

    for paired_align in tagged_paired_aligns:
        (left_umi, right_umi) =  paired_align.umi
        for best_tag in ranked_tags:
            if paired_align.umi == best_tag:
                tag_aligns[best_tag].add(paired_align)
                break
            elif left_umi == best_tag[0] or right_umi == best_tag[1]:
                tag_aligns[best_tag].add(paired_align)
                tag_inexact_match_count[best_tag] += 1
                break
            elif (_hamming_dist(left_umi, best_tag[0]) <= hamming_threshold) \
                or (_hamming_dist(right_umi, best_tag[1]) <= hamming_threshold):
                tag_aligns[best_tag].add(paired_align)
                tag_inexact_match_count[best_tag] += 1
                break
    tag_families = []
    for tag in sorted(tag_aligns):
        tag_family = _TagFamily(tag,
                               tag_aligns[tag],
                               tag_inexact_match_count[tag],
                               consensus_threshold,
                               min_family_size)
        tag_families.append(tag_family)
    #Necessary to make output deterministic
    tag_families.sort(key=lambda x: x.consensus.left_alignment.query_name)
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
    parser.add_argument("-d", "--umi_distance_threshold",
                        type=int,
                        default = DEFAULT_UMI_DISTANCE_THRESHOLD,
                        help=\
"""={} (>=0); UMIs equal to or closer than this Hamming distance will be
 combined into a single family. Lower threshold make more families with more
 consistent UMIs; 0 implies UMI must match
 exactly.""".format(DEFAULT_UMI_DISTANCE_THRESHOLD))

    args = parser.parse_args(arguments)
    return args

def _rank_tags(tagged_paired_aligns):
    '''Return the list of tags ranked from most to least popular.'''
    tag_count_dict = defaultdict(int)
    for paired_align in tagged_paired_aligns:
        umi =  paired_align.umi
        tag_count_dict[umi] += 1
    tags_by_count = sorted(tag_count_dict.items(),
                           key=lambda x: (-1 * x[1], x[0]))
    ranked_tags = [tag_count[0] for tag_count in tags_by_count]
    return ranked_tags


def _build_lightweight_pairs(aligned_segments, log):
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
    log.debug('filter_unpaired|{}/{} ({:.2f}%) alignments had mate present',
             len(lightweight_pairs) * 2,
             total_align_count,
             100 * len(lightweight_pairs) * 2 / total_align_count)
    log.debug('filter_unpaired|{}/{} ({:.2f}%) alignments missing their mate',
             len(name_pairs),
             total_align_count,
             100 * len(name_pairs) / total_align_count)
    return lightweight_pairs, total_align_count


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
            next_breakpoint = 10 * int(progress/10) + 10
        yield item
    log.info("100% ({}/{}) alignments processed", row_count, total_rows)

def _dedup_alignments(args, consensus_writer, annotated_writer, log):
    try:
        log.info('reading input bam [{}]', args.input_bam)
        bamfile = samtools.alignment_file(args.input_bam, 'rb')
        filtered_aligns = samtools.filter_alignments(bamfile.fetch(), log)
        (lightweight_pairs, 
         total_aligns) = _build_lightweight_pairs(filtered_aligns, log)
        bamfile.close()

        coord_manifest = _build_coordinate_read_name_manifest(lightweight_pairs)
        bamfile = samtools.alignment_file(args.input_bam, 'rb')

        handlers = familyhandler.build_family_handlers(args,
                                                       consensus_writer,
                                                       annotated_writer,
                                                       log)

        filtered_aligns_gen = samtools.filter_alignments(bamfile.fetch())
        progress_gen = _progress_logger(filtered_aligns_gen,
                                        total_aligns,
                                        log)
        for coord_family in _build_coordinate_families(progress_gen,
                                                       coord_manifest):
            ranked_tags = _rank_tags(coord_family)
            tag_families = _build_tag_families(coord_family,
                                               ranked_tags,
                                               args.umi_distance_threshold,
                                               args.consensus_freq_threshold,
                                               args.min_family_size_threshold)
            for handler in handlers:
                for tag_family in tag_families:
                    handler.handle(tag_family)

        for handler in handlers:
            handler.end()

        bamfile.close()

    except Exception: #pylint: disable=broad-except
        log.error("An unexpected error occurred")
        log.error(traceback.format_exc())
        exit(1)

def _log_environment_info(log, args):
    log.debug('original_command_line|{}',' '.join(args.original_command_line))
    log.debug('command_options|{}', vars(args))
    log.debug('command_cwd|{}', os.getcwd ())
    log.debug('platform_uname|{}', platform.uname())
    log.debug('platform_python_version|{}', platform.python_version())
    log.debug('pysam_version|{}', pysam.__version__)

def _build_writer(input_bam, output_bam, tags):
    if not output_bam:
        return samtools.AlignWriter.NULL
    else:
        input_bam = samtools.alignment_file(input_bam, "rb")
        header = input_bam.header
        input_bam.close()
        return samtools.AlignWriter(header, output_bam, tags)

def _build_bam_tags():
    tags = [
        samtools.BamTag("X0", "Z",
                        ("filter (rationale explaining why the align was "
                         "excluded)"),
                        lambda fam, align: align.filter),
        samtools.BamTag("X1", "i",
                        "unique identifier for this alignment family",
                        lambda fam, align: fam.umi_sequence if fam else None),
        samtools.BamTag("X2", "Z",
                        ("L~R UMT barcodes for this alignment family; because "
                         "of fuzzy matching the family UMT may be distinct "
                         "from the UMT of the original alignment"),
                        lambda fam, align: "{0}~{1}".format(fam.umi[0],
                                                            fam.umi[1]) if fam else None),
        samtools.BamTag("X3", "i",
                        "family size (number of align pairs in this family)",
                        lambda fam, align: len(fam.alignments) if fam else None),
        samtools.BamTag("X4", "i",
                        ("presence of this tag signals that this alignment "
                         "would be the template for the consensus alignment"),
                        lambda fam, align: 1 if fam and fam.consensus.left_alignment.query_name == align.query_name else None)]
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
        args = _parse_command_line_args(command_line_args[1:])
        args.original_command_line = command_line_args
        if not args.log_file:
            args.log_file = args.output_bam + ".log"
        log = utils.Logger(args)
        _log_environment_info(log, args)
        log.info('connor begins (v{})', __version__)
        log.info('logging to [{}]', args.log_file)
        bam_tags = _build_bam_tags()
        annotated_writer = _build_writer(args.input_bam,
                                         args.annotated_output_bam,
                                         bam_tags)
        consensus_writer = _build_writer(args.input_bam,
                                         args.output_bam,
                                         bam_tags)
        _dedup_alignments(args, consensus_writer, annotated_writer, log)
        annotated_writer.close()
        consensus_writer.close()
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


if __name__ == '__main__':
    main(sys.argv)
