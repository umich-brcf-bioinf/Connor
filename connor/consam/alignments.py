'''Connor-specific classes for working with alignments'''
from __future__ import print_function, absolute_import, division

import connor.utils as utils

DEFAULT_TAG_LENGTH = 6

class PairedAlignment(object):
    '''Represents the left and right aligns from a single paired sequence.'''
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

    def cigars(self, format_string=None):
        if format_string:
            return format_string.format(left=self.left.cigarstring,
                                        right=self.right.cigarstring)
        else:
            return self.left.cigarstring, self.right.cigarstring

    def positions(self, format_string=None):
        left_value = self.left.reference_start + 1
        right_value = self.right.reference_end + 1
        if format_string:
            return format_string.format(left=left_value, right=right_value)
        else:
            return left_value, right_value

    def replace_umt(self, umt):
        if not (umt[0] or umt[1]) or \
            (len(umt[0]) != self._tag_length) or \
            (len(umt[1]) != self._tag_length):
            msg = "Each UMT must match tag_length ({})"
            raise ValueError(msg.format(self._tag_length))
        left_qual = self.left.query_qualities
        right_qual = self.right.query_qualities
        left_query_frag = self.left.query_sequence[len(umt[0]):]
        left_query_frag_str = utils._byte_array_to_string(left_query_frag)
        self.left.query_sequence = umt[0] + left_query_frag_str
        right_query_frag = self.right.query_sequence[:-len(umt[1])]
        right_query_frag_str = utils._byte_array_to_string(right_query_frag)
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


class ConnorAlign(object):
    '''Wraps pysam alignment adding ability to track filtering criteria'''
    # cgates: FYI, you can use dynamic delegation via __setattr__ and
    # __getattr__ but it's awkward and about twice as slow
    def __init__(self, pysam_align_segment, filter_value=None):
        self.pysam_align_segment = pysam_align_segment
        self.filter_value = filter_value

    def __eq__(self, other):
        return other.__dict__ == self.__dict__

    # cgates: the native pysam hashing is not performant for ultradeep pileups
    def __hash__(self):
        return hash(self.filter_value) ^ \
               hash(self.pysam_align_segment.query_name) ^ \
               self.pysam_align_segment.reference_start

    @property
    def cigarstring(self):
        return self.pysam_align_segment.cigarstring
    @cigarstring.setter
    def cigarstring(self, value):
        self.pysam_align_segment.cigarstring = value

    @property
    def flag(self):
        return self.pysam_align_segment.flag
    @flag.setter
    def flag(self, value):
        self.pysam_align_segment.flag = value

    def get_tag(self, name, with_value_type=False):
        return self.pysam_align_segment.get_tag(name, with_value_type)

    def get_tags(self, with_value_type=False):
        return self.pysam_align_segment.get_tags(with_value_type)

    @property
    def mapping_quality(self):
        return self.pysam_align_segment.mapping_quality
    @mapping_quality.setter
    def mapping_quality(self, value):
        self.pysam_align_segment.mapping_quality = value

    @property
    def next_reference_start(self):
        return self.pysam_align_segment.next_reference_start
    @next_reference_start.setter
    def next_reference_start(self, value):
        self.pysam_align_segment.next_reference_start = value

    @property
    def orientation(self):
        if self.reference_start < self.next_reference_start:
            return 'left'
        elif self.reference_start > self.next_reference_start:
            return 'right'
        else:
            return 'neither'

    @property
    def query_name(self):
        return self.pysam_align_segment.query_name
    @query_name.setter
    def query_name(self, value):
        self.pysam_align_segment.query_name = value

    @property
    def query_sequence(self):
        return self.pysam_align_segment.query_sequence
    @query_sequence.setter
    def query_sequence(self, value):
        self.pysam_align_segment.query_sequence = value

    @property
    def query_qualities(self):
        return self.pysam_align_segment.query_qualities
    @query_qualities.setter
    def query_qualities(self, value):
        self.pysam_align_segment.query_qualities = value

    @property
    def reference_end(self):
        return self.pysam_align_segment.reference_end

    @property
    def reference_id(self):
        return self.pysam_align_segment.reference_id
    @reference_id.setter
    def reference_id(self, value):
        self.pysam_align_segment.reference_id = value

    @property
    def reference_name(self):
        return self.pysam_align_segment.reference_name

    @property
    def reference_start(self):
        return self.pysam_align_segment.reference_start
    @reference_start.setter
    def reference_start(self, value):
        self.pysam_align_segment.reference_start = value

    def set_tag(self, tag_name, tag_value, value_type):
        self.pysam_align_segment.set_tag(tag_name, tag_value, value_type)

    @property
    def template_length(self):
        return self.pysam_align_segment.template_length
    @template_length.setter
    def template_length(self, value):
        self.pysam_align_segment.template_length = value
