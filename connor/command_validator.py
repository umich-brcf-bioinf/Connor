'''Validates command preconditions.
Specifically checks that the command, arguments, and environment
(e.g. input/output directories or files) are consistent and plausible.
Each validation function evaluates a specific precondition raising a
UsageException if things look problematic or logging a warning if the condition
can be forced.
'''
from __future__ import print_function, absolute_import, division
import itertools
import os
import connor.utils as utils
import connor.samtools as samtools

_SAMPLE_SIZE = 1000

def _log_force_or_raise(args, log, msg):
    if args.force:
        log.warning(msg + ' (**forcing**)')
    else:
        raise utils.UsageError(msg + ' Are you sure? (--force to proceed)')

def _check_input_bam_exists(args, log=None): #pylint: disable=unused-argument
    if not os.path.exists(args.input_bam):
        raise utils.UsageError(("Specified input [{}] does not exist. Review "
                                "inputs and try again.").format(args.input_bam))

def _check_input_bam_valid(args, log=None): #pylint: disable=unused-argument
    try:
        bamfile = samtools.alignment_file(args.input_bam, 'rb')
        bamfile.close()
    except ValueError:
        raise utils.UsageError(("Specified input [{}] not a valid BAM. Review "
                                "inputs and try again.").format(args.input_bam))

def _check_input_bam_indexed(args, log=None): #pylint: disable=unused-argument
    bamfile = samtools.alignment_file(args.input_bam, 'rb')
    try:
        bamfile.fetch()
    except ValueError:
        raise utils.UsageError(("Specified input [{}] is not indexed. Review "
                                "inputs and try again.").format(args.input_bam))
    finally:
        bamfile.close()

def _check_input_bam_not_deduped(args, log=None):
    bamfile = samtools.alignment_file(args.input_bam, 'rb')
    header = bamfile.header
    bamfile.close()
    names = set([pg_item.get('PN', None) for pg_item in header.get('PG', [])])
    if samtools.CONNOR_PG_PN in names:
        msg = ('Specified input [{}] has already been processed with '
               'Connor.').format(args.input_bam)
        _log_force_or_raise(args, log, msg)

def _check_input_bam_not_empty(args, log=None): #pylint: disable=unused-argument
    bamfile = samtools.alignment_file(args.input_bam, 'rb')
    try:
        next(bamfile.fetch())
    except StopIteration:
        msg = "Specified input [{}] is empty"
        raise utils.UsageError(msg.format(args.input_bam))
    finally:
        bamfile.close()

def _check_input_bam_paired(args, log=None): #pylint: disable=unused-argument
    bamfile = samtools.alignment_file(args.input_bam, 'rb')
    try:
        for alignment in itertools.islice(bamfile.fetch(), _SAMPLE_SIZE):
            if alignment.is_paired:
                return
    finally:
        bamfile.close()
    msg = ('Specified input [{}] does not appear to contain paired '
           'reads').format(args.input_bam)
    _log_force_or_raise(args, log, msg)


def _check_input_bam_barcoded(args, log=None):
    #pylint: disable=invalid-name
    SOFTCLIP_OP = 4
    SOFTCLIP_THRESHOLD = 0.80

    def is_edge_softclipped(align):
        edge_index = -1 if align.is_reverse else 0
        return int(align.cigartuples[edge_index][0] == SOFTCLIP_OP)

    def strand(align):
        return 'reverse' if align.is_reverse else 'forward'

    softclipped = {'forward': [], 'reverse': []}
    bamfile = samtools.alignment_file(args.input_bam, 'rb')
    try:
        for align in itertools.islice(bamfile.fetch(), _SAMPLE_SIZE):
            softclipped[strand(align)].append(is_edge_softclipped(align))
    finally:
        bamfile.close()

    forward_freq = sum(softclipped['forward'])/len(softclipped['forward'])
    reverse_freq = sum(softclipped['reverse'])/len(softclipped['reverse'])
    if min(forward_freq, reverse_freq) < SOFTCLIP_THRESHOLD:
        msg = ('Specified input [{}] reads do not appear to have '
               'barcodes').format(args.input_bam)
        _log_force_or_raise(args, log, msg)

_VALIDATIONS = [_check_input_bam_exists,
                _check_input_bam_valid,
                _check_input_bam_indexed,
                _check_input_bam_not_deduped,
                _check_input_bam_not_empty,
                _check_input_bam_paired]


def preflight(args, log):
    for validate in _VALIDATIONS:
        validate(args, log)
        