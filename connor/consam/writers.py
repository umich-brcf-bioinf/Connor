'''Classes to simplify writing alignments to BAMs'''
from __future__ import print_function, absolute_import, division
from collections import OrderedDict
from collections import defaultdict
from copy import deepcopy

import connor
import connor.consam.pysamwrapper as pysamwrapper
import connor.utils as utils

class AlignWriter(object):
    '''Given header, and path, makes BAM; appends Connor headers/tags and indexes on close'''
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


CONNOR_PG_ID='connor'
CONNOR_PG_PN='connor'
_HEADER_PG_KEY = 'PG'

def _set_pg_header(header, simplify_pg_header, command_line):
    connor_header = {'ID':CONNOR_PG_ID,
                     'PN':CONNOR_PG_PN}
    if not simplify_pg_header:
        connor_header['CL'] = ' '.join(command_line)
        connor_header['VN'] = connor.__version__
    pg_headers = header.get(_HEADER_PG_KEY, [])
    pg_headers.append(connor_header)
    header[_HEADER_PG_KEY]=pg_headers

def build_writer(input_bam, output_bam, tags, args):
    if not output_bam:
        return AlignWriter.NULL

    input_bam = pysamwrapper.alignment_file(input_bam, 'rb')
    header_dict = pysamwrapper.get_header_dict(input_bam)
    input_bam.close()
    _set_pg_header(header_dict,
                   args.simplify_pg_header,
                   args.original_command_line)
    return AlignWriter(header_dict, output_bam, tags)
