#! /usr/bin/env python
'''Utils/classes for reading aligments'''
from __future__ import print_function, absolute_import, division
from collections import defaultdict

from connor.consam.bamflag import BamFlag
import connor.consam.writers as writers
import connor.consam.alignments as alignments
import connor.utils as utils

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
def _build_coordinate_pairs(connor_alignments, excluded_writer):
    MISSING_MATE_FILTER = 'read mate was missing or excluded'
    coords = defaultdict(dict)
    for alignment in connor_alignments:
        if alignment.orientation == 'left':
            key = (alignment.reference_id, alignment.next_reference_start)
            coords[key][alignment.query_name] = alignment
        elif alignment.orientation == 'neither':
            key = (alignment.reference_id, alignment.next_reference_start)
            if key in coords and alignment.query_name in coords[key]:
                align1 = coords[key].pop(alignment.query_name)
                yield alignments.PairedAlignment(align1, alignment)
            else:
                coords[key][alignment.query_name] = alignment
        else:
            key = (alignment.reference_id, alignment.reference_start)
            coord = coords[key]
            l_align = coord.pop(alignment.query_name, None)
            # Clear empty coordinate dict
            if not len(coord):
                del coords[key]
            if l_align:
                yield alignments.PairedAlignment(l_align, alignment)
            else:
                alignment.filter_value = MISSING_MATE_FILTER
                excluded_writer.write(None, None, alignment)
    for aligns in coords.values():
        for align in aligns.values():
            align.filter_value = MISSING_MATE_FILTER
            excluded_writer.write(None, None, align)

def _progress_logger(base_generator,
                     total_rows,
                     log,
                     supplemental_log=lambda x: None):
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
            supplemental_log(log)
            next_breakpoint = 10 * int(progress/10) + 10
        yield item
    log.info("100% ({}/{}) alignments processed", row_count, total_rows)
    supplemental_log(log)

def paired_reader(bamfile,
                  total_aligns,
                  log,
                  supplemental_log,
                  annotated_writer):
        progress_gen = _progress_logger(bamfile.fetch(),
                                        total_aligns,
                                        log,
                                        supplemental_log)
        filtered_aligns_gen = _filter_alignments(progress_gen,
                                                 annotated_writer)
        paired_align_gen = _build_coordinate_pairs(filtered_aligns_gen,
                                                   annotated_writer)
        return paired_align_gen
