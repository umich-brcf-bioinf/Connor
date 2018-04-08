'''Classes to help define tag and coordinate families.'''
from __future__ import print_function, absolute_import, division
try:
    from builtins import range as iter_range
except ImportError:
    from __builtin__ import xrange as iter_range
from collections import defaultdict
from collections import Counter
from copy import deepcopy
from functools import partial

from sortedcontainers import SortedSet

from connor.consam.alignments import PairedAlignment

class TagFamily(object):
    '''Holds alignments that share same coordinate and UMI (allowing fuzzy UMI matches)'''
    umi_sequence = 0

    def __init__(self,
                 umt,
                 alignments,
                 inexact_match_count,
                 consensus_threshold,
                 family_filter=lambda _: None):
        self.umi_sequence = TagFamily.umi_sequence
        TagFamily.umi_sequence += 1
        self._umt = umt
        (self.distinct_cigar_count,
         majority_cigar) = TagFamily._get_dominant_cigar_stats(alignments)
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
            if TagFamily._get_cigarstring_tuple(pair) != majority_cigar:
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
        counter = Counter([TagFamily._get_cigarstring_tuple(s) for s in alignments])
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
        template_pair = TagFamily._select_template_alignment_pair(included_pairs)

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
        consensus_pair = PairedAlignment(left_align,
                                                   right_align,
                                                   tag_length=len(umt[0]))
        consensus_pair.replace_umt(umt)
        return consensus_pair

class CoordinateFamilyHolder(object):
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
