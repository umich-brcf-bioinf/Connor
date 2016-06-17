"""Simplifies discrepancies in how different versions pysam wrap samtools"""
from __future__ import print_function, absolute_import, division

import pysam
import os
import re

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

def sort(input_bam_filepath, output_bam_filepath):
    SAMTOOLS_UTIL.sort(input_bam_filepath, output_bam_filepath)

def index(bam_filepath):
    SAMTOOLS_UTIL.index(bam_filepath)
