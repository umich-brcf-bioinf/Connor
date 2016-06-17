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

import os
import sys
import traceback
import pysam
import connor.samtools

__version__ = connor.__version__

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


DEFAULT_TAG_LENGTH = 6

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
        if pos1 < pos2:
            self.key = (chrom, pos1, pos2)
        else:
            self.key = (chrom, pos2, pos1)


class PairedAlignment(object):
    '''Represents the left and right align pairs from an single sequence.'''
    def __init__(self, left_alignment,
                 right_alignment,
                 tag_length=DEFAULT_TAG_LENGTH):
        self.left_alignment = left_alignment
        self.right_alignment = right_alignment
        self._tag_length = tag_length
        left_tag_id = self.left_alignment.query_sequence[0:self._tag_length]
        right_tag_id = self.right_alignment.query_sequence[0:self._tag_length]
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
    def __init__(self, umi, list_of_alignments):
        self.umi = umi
        self.alignments = list_of_alignments
        self.consensus = self._build_consensus(umi, list_of_alignments)

    @staticmethod
    def _generate_consensus_sequence(list_of_alignments):
        consensus = []
        for i in xrange(0, len(list_of_alignments[0].query_sequence)):
            counter = Counter([s.query_sequence[i:i+1] for s in list_of_alignments])
            base = counter.most_common(1)[0][0]
#            freq = 100 * counter[base]/sum(counter.values())
            consensus.append(base)
#             if freq >= consensus_cutoff:
#                 consensus.append(base)
#             else:
#                 consensus.append("N")
        return "".join(consensus)
    
    @staticmethod
    def _generate_consensus_qualities(list_of_alignments):
        consensus_quality = []
        for i in xrange(0, len(list_of_alignments[0].query_qualities)):
            counter = Counter([s.query_qualities[i:i+1] for s in list_of_alignments])
            base = counter.most_common(1)[0][0]
#            freq = 100 * counter[base]/sum(counter.values())
            consensus_quality.append(base)
#             if freq >= consensus_cutoff:
#                 consensus.append(base)
#             else:
#                 consensus.append("N")
        return "".join(consensus_quality)

    #TODO: (cgates) I don't like that the tags assume umi is a tuple and symmetric between left and right 
    @staticmethod
    def _build_consensus(umi, alignments):
        left_aligns = []
        right_aligns = []
        for align in alignments:
            left_aligns.append(align.left_alignment)
            right_aligns.append(align.right_alignment)
        left_consensus_sequence = TagFamily._generate_consensus_sequence(left_aligns)
        right_consensus_sequence = TagFamily._generate_consensus_sequence(right_aligns)
        left_consensus_qualities = TagFamily._generate_consensus_qualities(left_aligns)
        right_consensus_qualities = TagFamily._generate_consensus_qualities(right_aligns)
        left_consensus_align = deepcopy(alignments[0].left_alignment,
                                        {'query_sequence':left_consensus_sequence})
        right_consensus_align = deepcopy(alignments[0].right_alignment,
                                         {'query_sequence':right_consensus_sequence})
        left_consensus_align.query_qualities = left_consensus_qualities
        right_consensus_align.query_qualities = right_consensus_qualities
        consensus_align = PairedAlignment(left_consensus_align,
                                          right_consensus_align,
                                          tag_length=len(umi[0]))
        consensus_align.replace_umi(umi)
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
            key = LightweightAlignment(aseg).key
            family_dict[key].add(paired_align)
            coord_read_name_manifest[key].remove(aseg.query_name)
            if not coord_read_name_manifest[key]:
                yield family_dict.pop(key)

def _build_consensus_pair(alignment_family):
    '''Aggregate a set of reads into a single consensus read.'''
    sorted_reads = sorted(alignment_family,
                          key=lambda x: x.left_alignment.query_name)
    return sorted_reads.pop()


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
    #Sorting by tag is not strictly necessary but keeps results deterministic
    sorted_values = [tag_aligns[key] for key in sorted(tag_aligns)]
    return sorted_values

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
        lw_aligns = [LightweightAlignment(align) for align in bamfile.fetch()]
        original_read_count = len(lw_aligns)
        _log('original read count: {}', original_read_count)
        coord_manifest = _build_coordinate_read_name_manifest(lw_aligns)
        bamfile.close()
        bamfile = pysam.AlignmentFile(args.input_bam, 'rb')
        outfile = pysam.AlignmentFile(args.output_bam, 'wb', template=bamfile)
        consensus_read_count = 0
        for coord_family in _build_coordinate_families(bamfile.fetch(),
                                                       coord_manifest):
            ranked_tags = _rank_tags(coord_family)
            for tag_family in _build_tag_families(coord_family, ranked_tags):
                read_pair = _build_consensus_pair(tag_family)
                outfile.write(read_pair.left_alignment)
                outfile.write(read_pair.right_alignment)
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
