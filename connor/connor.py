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
from datetime import datetime
import itertools
import operator
import os
import sys
import traceback

import pandas as pd
import pysam
import connor.samtools

__version__ = connor.__version__

try:
    xrange
except NameError:
    xrange = range

try:
    iter_map = itertools.imap
except AttributeError:
    iter_map = map


DEFAULT_TAG_LENGTH = 6
DEFAULT_CONSENSUS_THRESHOLD=0.6
MIN_ORIG_READS = 3
HAMMING_THRESHOLD = 1


DESCRIPTION=\
'''Deduplicates BAM file based on custom inline DNA barcodes.
Emits a new BAM file reduced to a single consensus read for each family of
original reads.
'''

class _ConnorUsageError(Exception):
    """Raised for malformed command or invalid arguments."""
    def __init__(self, msg, *args):
        super(_ConnorUsageError, self).__init__(msg, *args)


class _ConnorArgumentParser(argparse.ArgumentParser):
    """Argument parser that raises UsageError instead of exiting."""
    #pylint: disable=too-few-public-methods
    def error(self, message):
        '''Suppress default exit behavior'''
        raise _ConnorUsageError(message)

#TODO: (cgates): Switch to Logger class w/ info,debug methods and verbose state
def _log(msg_format, *args):
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    try:
        print("{}|{}".format(timestamp, msg_format).format(*args),
              file=sys.stderr)
    except IndexError:
        print(args)
    sys.stderr.flush()


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
        left_start = min(aligned_segment1.reference_start, aligned_segment2.reference_start)
        right_end = max(aligned_segment1.reference_end, aligned_segment2.reference_end)
        self.key = (chrom, left_start, right_end)


class _PairedAlignment(object):
    '''Represents the left and right align pairs from an single sequence.'''
    def __init__(self, left_alignment,
                 right_alignment,
                 tag_length=DEFAULT_TAG_LENGTH):
        self.left_alignment = left_alignment
        self.right_alignment = right_alignment
        self._tag_length = tag_length
        left_tag_id = self.left_alignment.query_sequence[0:self._tag_length]
        right_tag_id = self.right_alignment.query_sequence[-1 * self._tag_length:]
        self.umi = (left_tag_id, right_tag_id)



    def replace_umi(self, umi):
        def byte_array_to_string(sequence):
            if isinstance(sequence, str):
                return sequence
            else:
                return str(sequence.decode("utf-8"))
        self.left_alignment.query_sequence = umi[0] + byte_array_to_string(self.left_alignment.query_sequence[len(umi[0]):])
        seq = byte_array_to_string(self.right_alignment.query_sequence)
        self.right_alignment.query_sequence = seq[:-1 * len(umi[1])] + umi[1]
        self.umi = umi

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


class _TagFamily(object):
    def __init__(self,
                 umi,
                 list_of_alignments,
                 inexact_match_count,
                 consensus_threshold=DEFAULT_CONSENSUS_THRESHOLD):
        self.umi = umi
        self.input_alignment_count = len(list_of_alignments)
        (self.distinct_cigar_count,
         self.minority_cigar_percentage,
         dominant_cigar) = _TagFamily._generate_dominant_cigar_stats(list_of_alignments)
        (self.alignments,
         self.excluded_alignments) = _TagFamily._get_alignments_for_dominant_cigar(dominant_cigar,
                                                                       list_of_alignments)
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
                excluded_alignments.append(pair)
        return included_alignments, excluded_alignments

    def _generate_consensus_sequence(self, list_of_alignments):
        consensus = []
        for i in xrange(0, len(list_of_alignments[0].query_sequence)):
            counter = Counter([s.query_sequence[i:i+1] for s in list_of_alignments])
            base = counter.most_common(1)[0][0]
            freq = counter[base] / sum(counter.values())
            if freq >= self.consensus_threshold:
                consensus.append(base)
            else:
                consensus.append("N")
        return "".join(consensus)

    #TODO (cgates): consider zipping over alignmnets instead
    @staticmethod
    def _generate_consensus_qualities(list_of_alignments):
        consensus_quality = []
        for i in xrange(0, len(list_of_alignments[0].query_qualities)):
            qualities = tuple([s.query_qualities[i] for s in list_of_alignments])
            counter = Counter(qualities)
            qual = counter.most_common(1)[0][0]
            consensus_quality.append(qual)
        return consensus_quality
    
    @staticmethod
    def _generate_dominant_cigar_stats(list_of_alignments):
        counter = Counter([_TagFamily._get_cigarstring_tuple(s) for s in list_of_alignments])
        number_distict_cigars = len(counter)
        top_two_cigar_count = counter.most_common(2)
        if len(top_two_cigar_count) == 1:
            minority_cigar_percentage = 0
        else:
            minority_cigar_percentage = top_two_cigar_count[1][1]/len(list_of_alignments)
        dominant_cigar = top_two_cigar_count[0][0]
        return number_distict_cigars, minority_cigar_percentage, dominant_cigar


    #TODO: (cgates) tags should not assume umi is a tuple and symmetric between left and right 
    def _build_consensus(self, umi, alignments):
        left_aligns = []
        right_aligns = []
        for align in alignments:
            left_aligns.append(align.left_alignment)
            right_aligns.append(align.right_alignment)
        left_consensus_sequence = self._generate_consensus_sequence(left_aligns)
        right_consensus_sequence = self._generate_consensus_sequence(right_aligns)
        left_consensus_qualities = _TagFamily._generate_consensus_qualities(left_aligns)
        right_consensus_qualities = _TagFamily._generate_consensus_qualities(right_aligns)
        left_consensus_align = deepcopy(alignments[0].left_alignment, {})
        right_consensus_align = deepcopy(alignments[0].right_alignment, {})
        left_consensus_align.query_sequence = left_consensus_sequence
        right_consensus_align.query_sequence = right_consensus_sequence
        consensus_align = _PairedAlignment(left_consensus_align,
                                          right_consensus_align,
                                          tag_length=len(umi[0]))
        consensus_align.replace_umi(umi)
        left_consensus_align.query_qualities = left_consensus_qualities
        right_consensus_align.query_qualities = right_consensus_qualities
        self._add_tags(consensus_align, len(alignments))
        return consensus_align

#TODO: (cgates): make this into a handler
    def _add_tags(self, consensus_align, num_alignments):
        x0 = "{0}|{1}".format(self.umi[0], self.umi[1])
        x2 = "{0},{1}".format(consensus_align.left_alignment.reference_start + 1,
                              consensus_align.right_alignment.reference_end)
        x3 = num_alignments
        x4 = num_alignments < MIN_ORIG_READS
        consensus_align.left_alignment.set_tag("X0", x0, "Z")
        consensus_align.left_alignment.set_tag("X2", x2, "Z")
        consensus_align.left_alignment.set_tag("X3", x3, "i")
        consensus_align.left_alignment.set_tag("X4", str(x4), "Z")

        consensus_align.right_alignment.set_tag("X0", x0, "Z")
        consensus_align.right_alignment.set_tag("X2", x2, "Z")
        consensus_align.right_alignment.set_tag("X3", x3, "i")
        consensus_align.right_alignment.set_tag("X4", str(x4), "Z")

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

def _build_coordinate_families(aligned_segments,coord_read_name_manifest):
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

def _build_tag_families(tagged_paired_aligns, ranked_tags):
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
            elif (_hamming_dist(left_umi, best_tag[0]) <= HAMMING_THRESHOLD) or \
                (_hamming_dist(right_umi, best_tag[1]) <= HAMMING_THRESHOLD):
                tag_aligns[best_tag].add(paired_align)
                tag_inexact_match_count[best_tag] += 1
                break
    tag_families = []
    for tag in tag_aligns:
        tag_family = _TagFamily(tag,
                               tag_aligns[tag],
                               tag_inexact_match_count[tag])
        tag_families.append(tag_family)
    #Necessary to make output deterministic
    tag_families.sort(key=lambda x: x.consensus.left_alignment.query_name)
    return tag_families

def _hamming_dist(str1, str2):
    assert len(str1) == len(str2)
    return sum(iter_map(operator.ne, str1, str2))

def _parse_command_line_args(arguments):
    parser = _ConnorArgumentParser( \
        formatter_class=argparse.RawTextHelpFormatter,
        usage="connor input_bam output_bam",
        description=(DESCRIPTION))

    parser.add_argument("-V",
                        "--version",
                        action='version',
                        version=__version__)
    parser.add_argument('input_bam',
                        help="path to input BAM")
    parser.add_argument('output_bam',
                        help="path to output BAM")
    args = parser.parse_args(arguments)
    return args

def _rank_tags(tagged_paired_aligns):
    '''Return the list of tags ranked from most to least popular.'''
    tag_count = defaultdict(int)
    for paired_align in tagged_paired_aligns:
        umi =  paired_align.umi
        tag_count[umi] += 1
    tags_by_count = sorted(tag_count.items(),
                           key=lambda x: (-1 * x[1], x[0]))
    ranked_tags = [tag_count[0] for tag_count in tags_by_count]
    return ranked_tags


def _sort_and_index_bam(bam_filename):
    output_dir = os.path.dirname(bam_filename)
    output_root = os.path.splitext(os.path.basename(bam_filename))[0]
    sorted_bam_filename = os.path.join(output_dir,
                                       output_root + ".sorted.bam")
    connor.samtools.sort(bam_filename, sorted_bam_filename)
    os.rename(sorted_bam_filename, bam_filename)
    connor.samtools.index(bam_filename)

def build_lightweight_aligns(aligned_segments):
    name_pairs = defaultdict(list)
    for align_segment in aligned_segments:
        name_pairs[align_segment.query_name].append(_LightweightAlignment(align_segment))
    return name_pairs
    
    
def build_lightweight_pairs(aligned_segments):
    name_pairs = dict()
    lightweight_pairs = list()
    for align_segment in aligned_segments:
        query_name = align_segment.query_name
        if not query_name in name_pairs: 
            name_pairs[align_segment.query_name] = align_segment
        else:
            new_pair = _LightweightPair(name_pairs.pop(query_name), align_segment)
            lightweight_pairs.append(new_pair)
    return lightweight_pairs

class _WriteFamilyHandler(object):
    def __init__(self, output_file,
                 output_filename,
                 min_original_pairs_theshold):
        self.output_file = output_file
        self.output_filename = output_filename
        self.min_original_pairs_theshold = min_original_pairs_theshold
        self.included_family_count = 0
        self.excluded_family_count = 0
        self.total_alignment_count = 0

    def handle(self, tag_families):
        for tag_family in tag_families:
            self.total_alignment_count += len(tag_family.alignments)
            if len(tag_family.alignments) >= self.min_original_pairs_theshold:
                consensus_pair = tag_family.consensus
                self.output_file.write(consensus_pair.left_alignment)
                self.output_file.write(consensus_pair.right_alignment)
                self.included_family_count += 1
            else:
                self.excluded_family_count += 1

    def end(self):
        total_family_count = self.excluded_family_count + \
                             self.included_family_count
        _log(('INFO|{}/{} ({:.2f}%) families were excluded because the '
              'original read count < {}'),
             self.excluded_family_count,
             total_family_count,
             100 * self.excluded_family_count / total_family_count,
             self.min_original_pairs_theshold)
        dedup_percent = 100 * (1 - (self.included_family_count / self.total_alignment_count))
        _log(('INFO|{} original pairs were deduplicated to {} families '
              '(overall dedup rate {:.2f}%)'),
             self.total_alignment_count,
             self.included_family_count,
             dedup_percent)
        _log('INFO|{} families written to [{}]',
             self.included_family_count,
             self.output_filename)


class _BaseTukeyStatHandler(object):
    def __init__(self):
        self.collection = []
        self.min = None
        self.quartile_1 = None
        self.median = None
        self.quartile_3 = None
        self.max = None

    def get_family_statistic(self, tag_family):
        pass

    def handle(self, tag_families):
        for tag_family in tag_families:
            stat = self.get_family_statistic(tag_family)
            if stat is not None:
                self.collection.append(stat)

    @property
    def summary(self):
        return (self.min,
                self.quartile_1,
                self.median,
                self.quartile_3,
                self.max)

    def end(self):
        summary = pd.Series(self.collection).describe()
        self.min = summary['min']
        self.max = summary['max']
        self.median = summary['50%']
        self.quartile_1 = summary['25%']
        self.quartile_3 = summary['75%']


class _FamilySizeStatHandler(_BaseTukeyStatHandler):
    def get_family_statistic(self, tag_family):
        return len(tag_family.alignments)

    def end(self):
        super(_FamilySizeStatHandler, self).end()
        _log('DEBUG|family stat|family size distribution (original pair counts: min, 1Q, median, 3Q, max): {}',
             ', '.join(map(str, self.summary)))


class _CigarMinorityStatHandler(_BaseTukeyStatHandler):
    def get_family_statistic(self, tag_family):
        stat = tag_family.minority_cigar_percentage
        if stat:
            return stat
        else:
            return None

    def end(self):
        super(_CigarMinorityStatHandler, self).end()
        _log('DEBUG|family stat|cigar|family distribution of minority CIGAR percentages (min, 1Q, median, 3Q, max): {}',
              ', '.join(map(lambda x: str(round(x,2)), self.summary)))
        
#Split into tukey and other stats
class _CigarStatHandler(object):
    def __init__(self):
        self.distinct_cigar_counts = []
        self.min = None
        self.quartile_1 = None
        self.median = None
        self.quartile_3 = None
        self.max = None
        self.total_input_alignment_count = 0
        self.total_alignment_count = 0
        self.total_family_count = 0
        self.families_cigar_counts = defaultdict(int)
        self.family_cigar_minority_percentage_counts = defaultdict(list)


    def handle(self, tag_families):
        for tag_family in tag_families:
            self.distinct_cigar_counts.append(tag_family.distinct_cigar_count)
            self.total_input_alignment_count += tag_family.input_alignment_count
            self.total_alignment_count += len(tag_family.alignments)
            self.total_family_count += 1
            self.families_cigar_counts[tag_family.distinct_cigar_count] += 1
            self.family_cigar_minority_percentage_counts[tag_family.distinct_cigar_count].append(tag_family.minority_cigar_percentage)


    @property
    def percent_deduplication(self):
        return 1 - (self.total_family_count / self.total_input_alignment_count)


    @property
    def total_excluded_alignments(self):
        return self.total_input_alignment_count - self.total_alignment_count

    @property
    def percent_excluded_alignments(self):
        return self.total_excluded_alignments / self.total_input_alignment_count

    @property
    def summary(self):
        return (self.min,
                self.quartile_1,
                self.median,
                self.quartile_3,
                self.max)

    def end(self):
        summary = pd.Series(self.distinct_cigar_counts).describe()
        self.min = summary['min']
        self.max = summary['max']
        self.median = summary['50%']
        self.quartile_1 = summary['25%']
        self.quartile_3 = summary['75%']

        _log('DEBUG|family stat|cigar|family distribution of distinct CIGAR counts (min, 1Q, median, 3Q, max): {}',
              ', '.join(map(str, self.summary)))
        ordered_cigar_counts = sorted(self.families_cigar_counts.items(),
                                      key = lambda x: -1 * int(x[1]))
        for num_cigars, freq in ordered_cigar_counts:
            summary = pd.Series(self.family_cigar_minority_percentage_counts[num_cigars]).describe()
            _log('DEBUG|family stat|cigar|{}/{} ({:.2f}%) families had {} CIGAR: minor % distrib {:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}',
                 freq,
                 self.total_family_count,
                 100 * freq / self.total_family_count,
                 num_cigars,
                 summary['min'],
                 summary['25%'],
                 summary['50%'],
                 summary['75%'],
                 summary['max'])
        _log('DEBUG|family stat|cigar|{}/{} ({:.2f}%) pairs were excluded as minority CIGAR',
             self.total_excluded_alignments,
             self.total_input_alignment_count,
             100 * self.percent_excluded_alignments)
        _log(('DEBUG|family stat|{} original pairs (of majority CIGAR) were deduplicated to {} families '
              '(majority CIGAR dedup rate {:.2f}%)'),
              self.total_alignment_count,
              self.total_family_count,
              100 * self.percent_deduplication)


class _MatchStatHandler(object):
    def __init__(self):
        self.total_inexact_match_count = 0
        self.total_pair_count = 0

    def handle(self, tag_families):
        for tag_family in tag_families:
            self.total_inexact_match_count += tag_family.inexact_match_count
            self.total_pair_count += len(tag_family.alignments)

    def end(self):
        _log(('DEBUG|family stat|{}/{} ({:.2f}%) original pairs matched by Hamming '
              'distance threshold (<={}) on left or right UMI '),
             self.total_inexact_match_count,
             self.total_pair_count,
             100 * self.percent_inexact_match,
             HAMMING_THRESHOLD)

    @property
    def percent_inexact_match(self):
        return self.total_inexact_match_count/self.total_pair_count

#TODO cgates: check that input file exists and output file does not
def main(command_line_args=None):
    '''Connor entry point.  See help for more info'''

    if not command_line_args:
        command_line_args = sys.argv

    try:
        args = _parse_command_line_args(command_line_args[1:])
        _log('INFO|connor begins (v{})', __version__)
        _log('INFO|reading input bam [{}]', args.input_bam)
        bamfile = pysam.AlignmentFile(args.input_bam, 'rb')
        lightweight_pairs = build_lightweight_pairs(bamfile.fetch())
        bamfile.close()

        original_read_count = len(lightweight_pairs) * 2
        _log('INFO|bam stat|{} original reads', original_read_count)
        coord_manifest = _build_coordinate_read_name_manifest(lightweight_pairs)
        bamfile = pysam.AlignmentFile(args.input_bam, 'rb')
        outfile = pysam.AlignmentFile(args.output_bam, 'wb', template=bamfile)

        handlers = [_FamilySizeStatHandler(),
                    _MatchStatHandler(),
                    _CigarMinorityStatHandler(),
                    _CigarStatHandler(),
                    _WriteFamilyHandler(outfile,
                                        args.output_bam,
                                        MIN_ORIG_READS)]

#TODO: (cgates): switch all handlers to family-at-at-time handling
        for coord_family in _build_coordinate_families(bamfile.fetch(),
                                                       coord_manifest):
            ranked_tags = _rank_tags(coord_family)
            tag_families = _build_tag_families(coord_family, ranked_tags)
            for handler in handlers:
                handler.handle(tag_families)

        for handler in handlers:
            handler.end()

        outfile.close()
        bamfile.close()

        _log('INFO|sorting and indexing bam')
        _sort_and_index_bam(args.output_bam)

        _log('INFO|connor complete')
    except _ConnorUsageError as usage_error:
        message = "connor usage problem: {}".format(str(usage_error))
        print(message, file=sys.stderr)
        print("See 'connor --help'.", file=sys.stderr)
        sys.exit(1)
    except Exception: #pylint: disable=broad-except
        _log("ERROR: An unexpected error occurred")
        _log(traceback.format_exc())
        exit(1)


if __name__ == '__main__':
    main(sys.argv)
