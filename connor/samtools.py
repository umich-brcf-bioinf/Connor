"""Simplifies discrepancies in how different versions pysam wrap samtools"""
from __future__ import print_function, absolute_import, division
from collections import defaultdict, OrderedDict
from copy import deepcopy
import os
import re
import pysam
import connor
import connor.utils as utils


class AlignWriter(object):
    class _NullWriter(object):
        def write(self, family, connor_align):
            pass

        def close(self, log=None):
            pass

    NULL= _NullWriter()
    BAM_TAG_FORMAT = 'connor|BAM tag|{}:{}'.replace('|', '\t')

    def __init__(self, header, bam_path, tags=None):
        if tags is None:
            self._tags = []
        else:
            self._tags = sorted(tags)
        new_header = self._add_header_lines(header, self._tags)
        self._bam_path = bam_path
        self._bam_file = pysam.AlignmentFile(bam_path, "wb", header=new_header)

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

    def _add_bam_tags(self, family, connor_align):
        for tag in self._tags:
            tag.set_tag(family, connor_align)

    def write(self, family, connor_align):
        self._add_bam_tags(family, connor_align)
        self._bam_file.write(connor_align.pysam_align_segment)

    def close(self, log=None):
        self._bam_file.close()
        if log:
            log.info("sorting/indexing [{}]".format(self._bam_path))
        sort_and_index_bam(self._bam_path)


class LoggingWriter(object):

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

    #TODO: cgates: promote UNPLACED_FAMILY to public NULL class; use throughout
    def write(self, family, connor_align):
        if not family:
            family = LoggingWriter.UNPLACED_FAMILY
        self._align_filter_stats[(family.filter_value,
                                  connor_align.filter_value)] += 1
        self._family_filter_stats[family.filter_value].add(family.umi_sequence)
        if family==LoggingWriter.UNPLACED_FAMILY:
            self._base_writer.write(None, connor_align)
        else:
            self._base_writer.write(family, connor_align)

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
        included_count = len(family_filter_stats.pop(None))
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

#     @staticmethod
#     def _log_stat(log_method, text, count, total):
#         log_method('{:.2f}% ({}/{}) {}',
#                    100 * count / total,
#                    count,
#                    total,
#                    text)

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
         total_align_count) = self._align_stats
        (included_fam_count,
         total_fam_count,
         discarded_fam_filter_counts) = self._family_stats

        self._log.info('{} alignments unplaced or discarded',
                       LoggingWriter._percent_stat_str(excluded_align_count,
                                                       total_align_count))

        LoggingWriter._log_filter_counts(self._unplaced_aligns,
                                         self._log.debug,
                                         'alignments unplaced: {percent_stat} {filter_name}',
                                         total_align_count)

        LoggingWriter._log_filter_counts(self._discarded_aligns,
                                         self._log.debug,
                                         'alignments discarded: {percent_stat} {filter_name}',
                                         total_align_count)

        LoggingWriter._log_filter_counts(discarded_fam_filter_counts,
                                         self._log.info,
                                         'families discarded: {percent_stat} {filter_name}',
                                         total_fam_count)

        percent_stat = LoggingWriter._percent_stat_str(included_align_count,
                                                       total_align_count)
        self._log.info('{} alignments included in {} families',
                       percent_stat,
                       included_fam_count)

        msg = ('{:.2f}% deduplication rate '
               '(1 - {} families/{} included alignments)')
        percent_dedup = 100 * (1 - (included_fam_count / included_align_count))
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

    def set_tag(self, family, connor_align):
        value = self._get_value(family, connor_align)
        connor_align.set_tag(self._tag_name, value, self._tag_type)


class BamFlag(object):
    PAIRED = 1
    PROPER_PAIR = 2
    UNMAP = 4
    MUNMAP = 8
    REVERSE = 16
    MREVERSE = 32
    READ1 = 64
    READ2 = 128
    SECONDARY = 256
    QCFAIL = 512
    DUP = 1024
    SUPPLEMENTARY = 2048

class _Pysam9SamtoolsUtil(object):
    @staticmethod
    def index(bam_filepath):
        pysam.samtools.index(bam_filepath, catch_stdout=False)

    @staticmethod
    def sort(input_bam_filepath, output_bam_filepath):
        pysam.samtools.sort(input_bam_filepath,
                            '-o',
                            output_bam_filepath,
                            catch_stdout=False)

    @staticmethod
    def _byte_array_to_string(sequence):
        if isinstance(sequence, str):
            return sequence
        else:
            return str(sequence.decode("utf-8"))

    @staticmethod
    def idxstats(input_bam_filepath):
        result = pysam.samtools.idxstats(input_bam_filepath)
        return _Pysam9SamtoolsUtil._byte_array_to_string(result).split('\n')


class _Pysam8SamtoolsUtil(object):
    @staticmethod
    def index(bam_filepath):
        pysam.index(bam_filepath, catch_stdout=False)

    @staticmethod
    def sort(input_bam_filepath, output_bam_filepath):
        output_bam_filepath_prefix = os.path.splitext(output_bam_filepath)[0]
        pysam.sort(input_bam_filepath,
                   output_bam_filepath_prefix,
                   catch_stdout=False)

    @staticmethod
    def idxstats(input_bam_filepath):
        return pysam.idxstats(input_bam_filepath)


def _get_samtools():
    if re.match(r".*\.9\.*", pysam.__version__):
        return _Pysam9SamtoolsUtil()
    else:
        return _Pysam8SamtoolsUtil()

SAMTOOLS_UTIL = _get_samtools()

class ConnorAlign(object):
    # cgates: FYI, you can use dynamic delegation via __setattr__ and
    # __getattr__ but it's awkward and about twice as slow
    def __init__(self, pysam_align_segment, filter_value=None):
        self.pysam_align_segment = pysam_align_segment
        self.filter_value = filter_value

    def __eq__(self, other):
        return other.__dict__ == self.__dict__

    def __hash__(self):
        return hash(self.filter_value) + hash(self.pysam_align_segment)

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
    filters = {'not in proper pair': \
                    lambda a: a.flag & BamFlag.PROPER_PAIR == 0,
                'secondary alignment': \
                    lambda a: a.flag & BamFlag.SECONDARY != 0,
                'qc failed': \
                    lambda a: a.flag & BamFlag.QCFAIL != 0,
                'mapping quality < 1': \
                    lambda a: a.mapping_quality < 1,
                'cigar unavailable': \
                    lambda a: a.cigarstring is None}
    generator = utils.FilteredGenerator(filters)
    for pysam_align, filter_value in generator.filter(pysam_alignments):
        connor_align = ConnorAlign(pysam_align, filter_value)
        if filter_value:
            excluded_writer.write(family=None,
                                  connor_align=connor_align)
        else:
            yield connor_align

def alignment_file(filename, mode, template=None):
    return pysam.AlignmentFile(filename, mode, template)

def sort(input_bam_filepath, output_bam_filepath):
    SAMTOOLS_UTIL.sort(input_bam_filepath, output_bam_filepath)

def index(bam_filepath):
    SAMTOOLS_UTIL.index(bam_filepath)

def sort_and_index_bam(bam_filename):
    output_dir = os.path.dirname(bam_filename)
    output_root = os.path.splitext(os.path.basename(bam_filename))[0]
    sorted_bam_filename = os.path.join(output_dir,
                                       output_root + ".sorted.bam")
    sort(bam_filename, sorted_bam_filename)
    os.rename(sorted_bam_filename, bam_filename)
    index(bam_filename)

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
        input_bam = alignment_file(input_bam, 'rb')
        header = input_bam.header
        input_bam.close()
        _set_pg_header(header,
                       args.simplify_pg_header,
                       args.original_command_line)
        return AlignWriter(header, output_bam, tags)

def total_align_count(input_bam):
    '''Returns count of all mapped alignments in input BAM (based on index)'''
    count = 0
    for line in SAMTOOLS_UTIL.idxstats(input_bam):
        if line:
            chrom, _, mapped, unmapped = line.strip().split('\t')
            if chrom != '*':
                count += int(mapped) + int(unmapped)
    return count
