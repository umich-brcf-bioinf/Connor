"""Simplifies discrepancies in how different versions pysam wrap samtools"""
from __future__ import print_function, absolute_import, division
from copy import deepcopy
import os
import re
import pysam
import connor.utils as utils


#TODO: cgates: make this into constants?
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


def _get_samtools():
    if re.match(r".*\.9\.*", pysam.__version__):
        return _Pysam9SamtoolsUtil()
    else:
        return _Pysam8SamtoolsUtil()

SAMTOOLS_UTIL = _get_samtools()

#TODO: cgates: This delegation pattern seems unclear
class ConnorAlign(object):
    def __init__(self, pysam_align_segment):
        self.__dict__['pysam_align_segment'] = pysam_align_segment
        self.__dict__['filter'] = None

    def __eq__(self, other):
        return other.__dict__ == self.__dict__

    def __hash__(self):
        return hash(self.__dict__['filter']) \
            + hash(self.__dict__['pysam_align_segment'])

    def __getattr__(self, name):
        return getattr(self.pysam_align_segment, name)

    def __setattr__(self, name, value):
        if name in self.__dict__:
            self.__dict__[name] = value
        else:
            delegator = self.__dict__['pysam_align_segment']
            delegator.__setattr__(name, value)

def filter_alignments(alignments, log=None):
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
    for alignment in generator.filter(alignments):
        yield ConnorAlign(alignment)
    if log:
        total = generator.total_excluded + generator.total_included
        log.debug(('filter_align|{}/{} ({:.2f}%) alignments passed filtering'),
                 generator.total_included,
                 total,
                 100 * generator.total_included / total)
        log.debug(('filter_align|{}/{} ({:.2f}%) alignments failed filtering '
                  'and will be excluded (see log file)'),
                 generator.total_excluded,
                 total,
                 100 * generator.total_excluded / total)
        for filter_name, count in generator.filter_stats.items():
            log.debug('filter_align|{:>7} alignments excluded because: {}',
                     count,
                     filter_name)

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


class AlignWriter(object):
    class _NullWriter(object):
        def write(self, *args):
            pass

        def close(self):
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

    def _add_bam_tags(self, family, connor_align):
        for tag in self._tags:
            tag.set_tag(family, connor_align)

    def write(self, family, connor_align):
        self._add_bam_tags(family, connor_align)
        self._bam_file.write(connor_align.pysam_align_segment)

    def close(self):
        self._bam_file.close()

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
