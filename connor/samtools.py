'''Connor-specific utils and classes for working with alignments'''
from __future__ import print_function, absolute_import, division
from collections import defaultdict, OrderedDict
from copy import deepcopy
import os

import connor
from connor.consam.bamflag import BamFlag
import connor.consam.pysamwrapper as pysamwrapper
import connor.utils as utils

DEFAULT_TAG_LENGTH = 6

class AlignWriter(object):
    class _NullWriter(object):
        def write(self, family, paired_align, connor_align):
            pass

        def close(self, log=None):
            pass

    NULL= _NullWriter()

    def __init__(self, header, bam_path, tags=None):
        if tags is None:
            self._tags = []
        else:
            self._tags = sorted(tags)
        new_header = self._add_header_lines(header, self._tags)
        self._bam_path = bam_path
        self._bam_file = pysamwrapper.alignment_file(bam_path,
                                                     'wb',
                                                     header=new_header)

    @staticmethod
    def _add_header_lines(original_header, tags):
        new_header =  deepcopy(original_header)
        if 'CO' not in new_header:
            new_header['CO'] = []
        new_header['CO'].extend([tag.header_comment for tag in tags])
        return new_header

    @property
    def bam_file_path(self):
        return self._bam_path

    def _add_bam_tags(self, family, paired_align, connor_align):
        for tag in self._tags:
            tag.set_tag(family, paired_align, connor_align)

    def write(self, family, paired_align, connor_align):
        self._add_bam_tags(family, paired_align, connor_align)
        self._bam_file.write(connor_align.pysam_align_segment)

    def close(self, log=None):
        self._bam_file.close()
        if log:
            log.info("sorting/indexing [{}]".format(self._bam_path))
        pysamwrapper.sort_and_index_bam(self._bam_path)


class LoggingWriter(object):
    #pylint: disable=too-few-public-methods
    class UnplacedFamily(object):
        def __init__(self):
            self.filter_value = 'unplaced'
            self.umi_sequence = -1

    UNPLACED_FAMILY = UnplacedFamily()

    def __init__(self, base_writer, log):
        self._base_writer = base_writer
        self._log = log
        self._align_filter_stats = defaultdict(int)
        self._family_filter_stats = defaultdict(set)

    def write(self, family, paired_align, connor_align):
        if not family:
            family = LoggingWriter.UNPLACED_FAMILY
        self._align_filter_stats[(family.filter_value,
                                  connor_align.filter_value)] += 1
        self._family_filter_stats[family.filter_value].add(family.umi_sequence)
        if family==LoggingWriter.UNPLACED_FAMILY:
            self._base_writer.write(None, paired_align, connor_align)
        else:
            self._base_writer.write(family, paired_align, connor_align)

    @staticmethod
    def _filter_counts(filter_dict):
        '''Returns an immutable ordered dict of filter:counts; when an item
        would be filtered by multiple filters, all are listed in alpha order;
        the dict itself is ordered by descending count, filter name.
        '''
        return OrderedDict(utils.sort_dict(filter_dict))

    @staticmethod
    def _log_line(text, count, total, filter_name):
        line = '{:.2f}% ({}/{}) {}: {}'
        return line.format(100 * count / total,
                           count,
                           total,
                           text,
                           filter_name)

    @property
    def _unplaced_aligns(self):
        unplaced_aligns = {}
        for (fam_filter, align_filter), cnt in self._align_filter_stats.items():
            if fam_filter == LoggingWriter.UNPLACED_FAMILY.filter_value and align_filter:
                unplaced_aligns[align_filter] = cnt
        return LoggingWriter._filter_counts(unplaced_aligns)

    @staticmethod
    def _discarded_filter_value(fam_filter, align_filter):
        if fam_filter == LoggingWriter.UNPLACED_FAMILY.filter_value:
            return None
        else:
            filter_values = []
            if fam_filter:
                filter_values.append(fam_filter)
            if align_filter:
                filter_values.append(align_filter)
            return "; ".join(filter_values)

    @property
    def _discarded_aligns(self):
        discarded_aligns = {}
        for (fam_filter, align_filter), cnt in self._align_filter_stats.items():
            filter_value = LoggingWriter._discarded_filter_value(fam_filter,
                                                                 align_filter)
            if filter_value:
                discarded_aligns[filter_value] = cnt
        return LoggingWriter._filter_counts(discarded_aligns)

    @property
    def _family_stats(self):
        family_filter_stats = dict(self._family_filter_stats)
        family_filter_stats.pop(LoggingWriter.UNPLACED_FAMILY.filter_value,
                                None)
        included_count = len(family_filter_stats.pop(None, []))
        discarded_count = 0
        filter_counts = OrderedDict()
        for name, fam_ids in family_filter_stats.items():
            align_count = len(fam_ids)
            discarded_count += align_count
            filter_counts[name] = align_count
        total_count = included_count + discarded_count
        return (included_count,
                total_count,
                LoggingWriter._filter_counts(filter_counts))

    @property
    def _align_stats(self):
        included_filter = (None, None)
        included_count = self._align_filter_stats[included_filter]
        excluded_count = sum([count for fam_align_filter, count in self._align_filter_stats.items() if fam_align_filter != included_filter])
        total_count = included_count + excluded_count
        return included_count, excluded_count, total_count

    @staticmethod
    def _percent_stat_str(count, total):
        return '{:.2f}% ({}/{})'.format(100 * count / total, count, total)

    @staticmethod
    def _log_filter_counts(filter_counts, log_method, msg_format, total):
        for name, count in filter_counts.items():
            percent = LoggingWriter._percent_stat_str(count, total)
            log_method(msg_format.format(filter_name=name,
                                         percent_stat=percent))

    def _log_results(self):
        (included_align_count,
         excluded_align_count,
         tot_align_count) = self._align_stats
        (included_fam_count,
         total_fam_count,
         discarded_fam_filter_counts) = self._family_stats

        self._log.info('{} alignments unplaced or discarded',
                       LoggingWriter._percent_stat_str(excluded_align_count,
                                                       tot_align_count))

        LoggingWriter._log_filter_counts(self._unplaced_aligns,
                                         self._log.debug,
                                         ('alignments unplaced: '
                                          '{percent_stat} {filter_name}'),
                                         tot_align_count)

        LoggingWriter._log_filter_counts(self._discarded_aligns,
                                         self._log.debug,
                                         ('alignments discarded: '
                                          '{percent_stat} {filter_name}'),
                                         tot_align_count)

        LoggingWriter._log_filter_counts(discarded_fam_filter_counts,
                                         self._log.info,
                                         ('families discarded: '
                                          '{percent_stat} {filter_name}'),
                                         total_fam_count)

        percent_stat = LoggingWriter._percent_stat_str(included_align_count,
                                                       tot_align_count)
        self._log.info('{} alignments included in {} families',
                       percent_stat,
                       included_fam_count)

        if included_align_count == 0:
            self._log.warning("No alignments passed filters. (Was input BAM downsampled?)")
        else:
            percent_dedup = 100 * (1 - (included_fam_count / included_align_count))
            msg = ('{:.2f}% deduplication rate '
                   '(1 - {} families/{} included alignments)')
            self._log.info(msg,
                           percent_dedup,
                           included_fam_count,
                           included_align_count)

    def close(self, log=None):
        if self._align_filter_stats:
            self._log_results()
        self._base_writer.close(log)


class BamTag(object):
    HEADER_FORMAT = 'connor|BAM tag|{}: {}'.replace('|', '\t')

    class _NullObject(object):
        '''Returns None for all method calls'''
        #pylint: disable=no-self-use, unused-argument
        def __init__(self):
            self.included_pair_count = None
            self.filter_value = None
            self.umi_sequence = None
            self.umt=lambda *args: None
            self.is_consensus_template = lambda *args: None
            self.positions = lambda *args: None
            self.cigars = lambda *args: None

    _NULL_OBJECT = _NullObject()
    def __init__(self, tag_name, tag_type, description, get_value):
        self._tag_name = tag_name
        self._tag_type = tag_type
        self._get_value = get_value
        self._description = description
        self.header_comment = BamTag.HEADER_FORMAT.format(tag_name,
                                                          description)

    def __lt__(self, other):
        return (self._tag_name,
                self._description) < (other._tag_name, other._description)

    def set_tag(self, family, paired_align, connor_align):
        family = family if family else BamTag._NULL_OBJECT
        paired_align = paired_align if paired_align else BamTag._NULL_OBJECT
        value = self._get_value(family, paired_align, connor_align)
        connor_align.set_tag(self._tag_name, value, self._tag_type)


class PairedAlignment(object):
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


def filter_alignments(pysam_alignments, excluded_writer=AlignWriter.NULL):
    filters = {'cigar unavailable': \
                    lambda a: a.cigarstring is None,
                'mapping quality < 1': \
                    lambda a: a.mapping_quality < 1,
                'not in proper pair': \
                    lambda a: a.flag & BamFlag.PROPER_PAIR == 0,
                'qc failed': \
                    lambda a: a.flag & BamFlag.QCFAIL != 0,
                'secondary alignment': \
                    lambda a: a.flag & BamFlag.SECONDARY != 0,
                'supplementary alignment': \
                    lambda a: a.flag & BamFlag.SUPPLEMENTARY != 0,
                    }
    generator = utils.FilteredGenerator(filters)
    for pysam_align, filter_value in generator.filter(pysam_alignments):
        connor_align = ConnorAlign(pysam_align, filter_value)
        if filter_value:
            excluded_writer.write(family=None,
                                  paired_align=None,
                                  connor_align=connor_align)
        else:
            yield connor_align


CONNOR_PG_ID='connor'
CONNOR_PG_PN='connor'
HEADER_PG_KEY = 'PG'

def _set_pg_header(header, simplify_pg_header, command_line):
    connor_header = {'ID':CONNOR_PG_ID,
                     'PN':CONNOR_PG_PN}
    if not simplify_pg_header:
        connor_header['CL'] = ' '.join(command_line)
        connor_header['VN'] = connor.__version__
    pg_headers = header.get(HEADER_PG_KEY, [])
    pg_headers.append(connor_header)
    header[HEADER_PG_KEY]=pg_headers

def build_writer(input_bam, output_bam, tags, args):
    if not output_bam:
        return AlignWriter.NULL
    else:
        input_bam = pysamwrapper.alignment_file(input_bam, 'rb')
        header_dict = pysamwrapper.get_header_dict(input_bam)
        input_bam.close()
        _set_pg_header(header_dict,
                       args.simplify_pg_header,
                       args.original_command_line)
        return AlignWriter(header_dict, output_bam, tags)

def total_align_count(input_bam):
    '''Returns count of all mapped alignments in input BAM (based on index)'''
    count = 0
    for line in pysamwrapper.idxstats(input_bam):
        if line:
            chrom, _, mapped, unmapped = line.strip().split('\t')
            if chrom != '*':
                count += int(mapped) + int(unmapped)
    return count

def _build_bam_tags():
    #pylint: disable=unused-argument
    def combine_filters(fam, paired_align, align):
        filters = [x.filter_value for x in [fam, align] if x and x.filter_value]
        if filters:
            return ";".join(filters).replace('; ', ';')
        else:
            return None
    boolean_tag_value = {True:1}
    tags = [
        BamTag("X0", "Z",
               ("filter (why the alignment was excluded)"),
               combine_filters),
        BamTag("X1", "Z",
               ("leftmost~rightmost matched pair positions"),
               lambda fam, pair, align: pair.positions('{left}~{right}')),
        BamTag("X2", "Z",
               ("L~R CIGARs"),
               lambda fam, pair, align: pair.cigars('{left}~{right}')),
        BamTag("X3", "i",
               "unique identifier for this alignment family",
               lambda fam, pair, align: fam.umi_sequence),
        BamTag("X4", "Z",
               ("L~R UMT barcodes for this alignment family; because "
                "of fuzzy matching the family UMT may be distinct "
                "from the UMT of the original alignment"),
               lambda fam, pair, align: fam.umt('{left}~{right}')),
        BamTag("X5", "i",
               "family size (number of align pairs in this family)",
               lambda fam, pair, align: fam.included_pair_count),
        BamTag("X6", "i",
               ("presence of this tag signals that this alignment "
                "would be the template for the consensus alignment"),
               lambda fam, pair, align: boolean_tag_value.get(fam.is_consensus_template(align), None))]
    return tags
