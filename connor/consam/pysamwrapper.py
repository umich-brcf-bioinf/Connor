"""Simplifies discrepancies in how different versions pysam wrap samtools"""
from __future__ import print_function, absolute_import, division
import os
from pkg_resources import parse_version

import pysam

import connor.utils as utils

class _Pysam8Wrapper(object):
    @staticmethod
    def is_compatible_version(pysam_version):
        return parse_version(pysam_version) < parse_version('0.9')

    @staticmethod
    def aligned_segment():
        return pysam.AlignedSegment()

    @staticmethod
    def get_header_dict(input_bam):
        return input_bam.header

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


class _Pysam9Wrapper(object):
    @staticmethod
    def is_compatible_version(pysam_version):
        return parse_version('0.9') <= parse_version(pysam_version) < parse_version('0.10')

    @staticmethod
    def aligned_segment():
        return pysam.AlignedSegment()

    @staticmethod
    def get_header_dict(input_bam):
        return input_bam.header

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
    def idxstats(input_bam_filepath):
        result = pysam.samtools.idxstats(input_bam_filepath)
        return utils._byte_array_to_string(result).split('\n')


class _Pysam10_11_12_13Wrapper(object):
    @staticmethod
    def is_compatible_version(pysam_version):
        return parse_version('0.10') <= parse_version(pysam_version) < parse_version('0.14')

    @staticmethod
    def aligned_segment():
        return pysam.AlignedSegment()

    @staticmethod
    def get_header_dict(input_bam):
        return input_bam.header

    @staticmethod
    def index(bam_filepath):
        pysam.samtools.index(bam_filepath, catch_stdout=False)

    @staticmethod
    def sort(input_bam_filepath, output_bam_filepath):
        pysam.samtools.sort('-o',
                            output_bam_filepath,
                            input_bam_filepath,
                            catch_stdout=False)

    @staticmethod
    def idxstats(input_bam_filepath):
        result = pysam.samtools.idxstats(input_bam_filepath)
        return utils._byte_array_to_string(result).split('\n')


class _Pysam14Wrapper(object):
    @staticmethod
    def is_compatible_version(pysam_version):
        return parse_version(pysam_version) >= parse_version('0.14')

    @staticmethod
    def aligned_segment():
        return pysam.AlignedSegment()

    @staticmethod
    def get_header_dict(input_bam):
        return input_bam.header.to_dict()

    @staticmethod
    def index(bam_filepath):
        pysam.samtools.index(bam_filepath, catch_stdout=False)

    @staticmethod
    def sort(input_bam_filepath, output_bam_filepath):
        pysam.samtools.sort('-o',
                            output_bam_filepath,
                            input_bam_filepath,
                            catch_stdout=False)

    @staticmethod
    def idxstats(input_bam_filepath):
        result = pysam.samtools.idxstats(input_bam_filepath)
        return utils._byte_array_to_string(result).split('\n')


def _get_pysam_wrapper():
    pysam_wrappers = [_Pysam14Wrapper(),
                      _Pysam10_11_12_13Wrapper(),
                      _Pysam9Wrapper(),
                      _Pysam8Wrapper()]
    for wrapper in pysam_wrappers:
        if wrapper.is_compatible_version(pysam.__version__):
            return wrapper
    msg = 'no wrapper compatible with pysam version {}'.format(pysam.__version__)
    raise RuntimeError(msg)

_WRAPPER = _get_pysam_wrapper()

def aligned_segment():
    return _WRAPPER.aligned_segment()

def alignment_file(filename, mode, header=None, template=None):
    if header:
        return pysam.AlignmentFile(filename, mode, header=header)
    if template:
        return pysam.AlignmentFile(filename, mode, template=template)
    return pysam.AlignmentFile(filename, mode)

def get_header_dict(bam):
    return _WRAPPER.get_header_dict(bam)

def idxstats(bam_filepath):
    return _WRAPPER.idxstats(bam_filepath)

def index(bam_filepath):
    _WRAPPER.index(bam_filepath)

def sort(input_bam_filepath, output_bam_filepath):
    _WRAPPER.sort(input_bam_filepath, output_bam_filepath)

def sort_and_index_bam(bam_filepath):
    output_dir = os.path.dirname(bam_filepath)
    output_root = os.path.splitext(os.path.basename(bam_filepath))[0]
    sorted_bam_filename = os.path.join(output_dir,
                                       output_root + ".sorted.bam")
    sort(bam_filepath, sorted_bam_filename)
    os.rename(sorted_bam_filename, bam_filepath)
    index(bam_filepath)

def total_align_count(bam_filepath):
    '''Returns count of all mapped alignments in input BAM (based on index)'''
    count = 0
    for line in idxstats(bam_filepath):
        if line:
            chrom, _, mapped, unmapped = line.strip().split('\t')
            if chrom != '*':
                count += int(mapped) + int(unmapped)
    return count
