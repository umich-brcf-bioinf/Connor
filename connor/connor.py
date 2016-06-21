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
import array
from collections import defaultdict, Counter
from copy import deepcopy
from datetime import datetime
import os
import sys
import traceback
import pysam
import connor.samtools

__version__ = connor.__version__

try:
    xrange
except NameError:
    xrange = range

DEFAULT_TAG_LENGTH = 6
DEFAULT_CONSENSUS_THRESHOLD=0.6
MIN_ORIG_READS = 3


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


def _log(msg_format, *args):
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    try:
        print("{}|{}".format(timestamp, msg_format).format(*args),
              file=sys.stderr)
    except IndexError:
        print(args)
    sys.stderr.flush()


class LightweightAlignment(object):
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


class LightweightPair(object):
    '''Minimal info from PySam.AlignedSegment used to expedite pos grouping.'''
    def __init__(self, aligned_segment1, aligned_segment2):
        self.name = aligned_segment1.query_name
        chrom = aligned_segment1.reference_name
        left_start = min(aligned_segment1.reference_start, aligned_segment2.reference_start)
        right_end = max(aligned_segment1.reference_end, aligned_segment2.reference_end)
        self.key = (chrom, left_start, right_end)


class PairedAlignment(object):
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
        self.left_alignment.query_sequence = umi[0] + self.left_alignment.query_sequence[len(umi[0]):]
        self.right_alignment.query_sequence = umi[1] + self.right_alignment.query_sequence[len(umi[1]):]
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


class TagFamily(object):
    def __init__(self,
                 umi,
                 list_of_alignments,
                 consensus_threshold=DEFAULT_CONSENSUS_THRESHOLD):
        self.umi = umi
        dominant_cigar = TagFamily._generate_dominant_cigar_pair(list_of_alignments)
        self.alignments = TagFamily._get_alignments_for_dominant_cigar(dominant_cigar,
                                                                       list_of_alignments)
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
        return [a for a in list_of_alignments if TagFamily._get_cigarstring_tuple(a) == dominant_cigar]

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
    def _generate_dominant_cigar_pair(list_of_alignments):
        counter = Counter([TagFamily._get_cigarstring_tuple(s) for s in list_of_alignments])
        return counter.most_common(1)[0][0]


    #TODO: (cgates) tags should not assume umi is a tuple and symmetric between left and right 
    def _build_consensus(self, umi, alignments):
        left_aligns = []
        right_aligns = []
        for align in alignments:
            left_aligns.append(align.left_alignment)
            right_aligns.append(align.right_alignment)
        left_consensus_sequence = self._generate_consensus_sequence(left_aligns)
        right_consensus_sequence = self._generate_consensus_sequence(right_aligns)
        left_consensus_qualities = TagFamily._generate_consensus_qualities(left_aligns)
        right_consensus_qualities = TagFamily._generate_consensus_qualities(right_aligns)
        left_consensus_align = deepcopy(alignments[0].left_alignment, {})
        right_consensus_align = deepcopy(alignments[0].right_alignment, {})
        left_consensus_align.query_sequence = left_consensus_sequence
        right_consensus_align.query_sequence = right_consensus_sequence
        consensus_align = PairedAlignment(left_consensus_align,
                                          right_consensus_align,
                                          tag_length=len(umi[0]))
        consensus_align.replace_umi(umi)
        left_consensus_align.query_qualities = left_consensus_qualities
        right_consensus_align.query_qualities = right_consensus_qualities
        self._add_tags(consensus_align, len(alignments))
        return consensus_align


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
            paired_align = PairedAlignment(pairing_dict.pop(aseg.query_name),
                                           aseg)
            key = LightweightPair(paired_align.left_alignment, 
                                  paired_align.right_alignment).key
            family_dict[key].add(paired_align)
            coord_read_name_manifest[key].remove(aseg.query_name)
            if not coord_read_name_manifest[key]:
                yield family_dict.pop(key)

def _build_tag_families(tagged_paired_aligns, ranked_tags):
    '''Return a list of read families; each family is a set of original reads.

    Each read is considered against each ranked tag until all reads are
    partitioned into families.'''
    tag_aligns = defaultdict(set)
    for paired_align in tagged_paired_aligns:
        (left_umi, right_umi) =  paired_align.umi
        for best_tag in ranked_tags:
            if left_umi == best_tag[0] or right_umi == best_tag[1]:
                tag_aligns[best_tag].add(paired_align)
                break
    tag_families = [TagFamily(tag, aligns) for tag, aligns in tag_aligns.items()]
    #Necessary to make output deterministic
    tag_families.sort(key=lambda x: x.consensus.left_alignment.query_name)
    return tag_families

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
        name_pairs[align_segment.query_name].append(LightweightAlignment(align_segment))
    return name_pairs
    
    
def build_lightweight_pairs(aligned_segments):
    name_pairs = dict()
    lightweight_pairs = list()
    for align_segment in aligned_segments:
        query_name = align_segment.query_name
        if not query_name in name_pairs: 
            name_pairs[align_segment.query_name] = align_segment
        else:
            new_pair = LightweightPair(name_pairs.pop(query_name), align_segment)
            lightweight_pairs.append(new_pair)
    return lightweight_pairs
    

#TODO cgates: check that input file exists and output file does not
def main(command_line_args=None):
    '''Connor entry point.  See help for more info'''

    if not command_line_args:
        command_line_args = sys.argv

    try:
        args = _parse_command_line_args(command_line_args[1:])
        _log('connor begins')
        _log('reading input bam  [{}]', args.input_bam)
        bamfile = pysam.AlignmentFile(args.input_bam, 'rb')
        lightweight_pairs = build_lightweight_pairs(bamfile.fetch())
        bamfile.close()

        original_read_count = len(lightweight_pairs) * 2
        _log('original read count: {}', original_read_count)
        coord_manifest = _build_coordinate_read_name_manifest(lightweight_pairs)
        bamfile = pysam.AlignmentFile(args.input_bam, 'rb')
        outfile = pysam.AlignmentFile(args.output_bam, 'wb', template=bamfile)
        consensus_read_count = 0
        for coord_family in _build_coordinate_families(bamfile.fetch(),
                                                       coord_manifest):
            ranked_tags = _rank_tags(coord_family)
            for tag_family in _build_tag_families(coord_family, ranked_tags):
                if len(tag_family.alignments) >= MIN_ORIG_READS:
                    consensus_pair = tag_family.consensus
                    outfile.write(consensus_pair.left_alignment)
                    outfile.write(consensus_pair.right_alignment)
                    consensus_read_count += 2
        _log('consensus read count: {}', consensus_read_count)
        _log('consensus/original: {:.4f}',
             consensus_read_count / original_read_count)
        outfile.close()
        bamfile.close()
    
        _log('sorting and indexing bam')
        _sort_and_index_bam(args.output_bam)
    
        _log('wrote deduped bam [{}]', args.output_bam)
        _log('connor complete')
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
