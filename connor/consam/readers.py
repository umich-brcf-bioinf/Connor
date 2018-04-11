'''Utils for reading aligments'''
from __future__ import print_function, absolute_import, division
from collections import defaultdict

from connor.consam.bamflag import BamFlag
import connor.consam.writers as writers
import connor.consam.alignments as alignments
import connor.utils as utils
from connor.consam import pysamwrapper

_MISSING_MATE_FILTER_TEXT = 'read mate was missing or excluded'

def _filter_alignments(pysam_alignments,
                      excluded_writer=writers.AlignWriter.NULL):
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
        connor_align = alignments.ConnorAlign(pysam_align, filter_value)
        if filter_value:
            excluded_writer.write(family=None,
                                  paired_align=None,
                                  connor_align=connor_align)
        else:
            yield connor_align

#TODO: cgates: reduce complexity
def _build_coordinate_pairs(umi_length, connor_alignments, excluded_writer):
    coords = defaultdict(dict)
    for alignment in connor_alignments:
        if alignment.orientation == 'left':
            key = (alignment.reference_id, alignment.next_reference_start)
            coords[key][alignment.query_name] = alignment
        elif alignment.orientation == 'neither':
            key = (alignment.reference_id, alignment.next_reference_start)
            if key in coords and alignment.query_name in coords[key]:
                align1 = coords[key].pop(alignment.query_name)
                yield alignments.PairedAlignment(align1, alignment, umi_length)
            else:
                coords[key][alignment.query_name] = alignment
        else:
            key = (alignment.reference_id, alignment.reference_start)
            coord = coords[key]
            l_align = coord.pop(alignment.query_name, None)
            # Clear empty coordinate dict
            if not coord:
                del coords[key]
            if l_align:
                yield alignments.PairedAlignment(l_align, alignment, umi_length)
            else:
                alignment.filter_value = _MISSING_MATE_FILTER_TEXT
                excluded_writer.write(None, None, alignment)
    for aligns in coords.values():
        for align in aligns.values():
            align.filter_value = _MISSING_MATE_FILTER_TEXT
            excluded_writer.write(None, None, align)

def _progress_logger(base_generator,
                     total_rows,
                     log,
                     log_usage=lambda: None):
    row_count = 0
    next_breakpoint = 0
    for item in base_generator:
        row_count += 1
        progress = 100 * row_count / total_rows
        if progress >= next_breakpoint and progress < 100:
            log.info("{}% ({}/{}) alignments processed",
                     next_breakpoint,
                     row_count,
                     total_rows)
            log_usage()
            next_breakpoint = 10 * int(progress/10) + 10
        yield item
    log.info("100% ({}/{}) alignments processed", row_count, total_rows)
    log_usage()

def _paired_reader(umi_length,
                   bamfile_gen,
                   total_aligns,
                   log,
                   supplemental_log,
                   annotated_writer):
    progress_gen = _progress_logger(bamfile_gen,
                                    total_aligns,
                                    log,
                                    supplemental_log)
    filtered_aligns_gen = _filter_alignments(progress_gen,
                                             annotated_writer)
    paired_align_gen = _build_coordinate_pairs(umi_length,
                                               filtered_aligns_gen,
                                               annotated_writer)
    return paired_align_gen

def _bamfile_generator(bam_filename):
    bamfile = pysamwrapper.alignment_file(bam_filename, 'rb')
    for align in bamfile.fetch():
        yield align
    bamfile.close()

def paired_reader_from_bamfile(args,
                               log,
                               usage_logger,
                               annotated_writer):
    '''Given a BAM file, return a generator that yields filtered, paired reads'''
    total_aligns = pysamwrapper.total_align_count(args.input_bam)
    bamfile_generator  = _bamfile_generator(args.input_bam)
    return _paired_reader(args.umt_length,
                          bamfile_generator,
                          total_aligns,
                          log,
                          usage_logger,
                          annotated_writer)
