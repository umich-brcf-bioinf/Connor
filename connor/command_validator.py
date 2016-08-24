'''Validates command preconditions.
Specifically checks that the command, arguments, and environment
(e.g. input/output directories or files) are consistent and plausible.
Each validation function evaluates a specific precondition.
Each function is allowed to:
 * change the environment (e.g. create a dir)
 * change the args (replace the original output dir with a new temp output dir)
 * raise a UsageException if things look problematic
'''
from __future__ import print_function, absolute_import, division

import os
import connor.utils as utils
import connor.samtools as samtools
import itertools

def _check_input_bam_exists(args, log=None):
    if not os.path.exists(args.input_bam):
        raise utils.UsageError(("Specified input [{}] does not exist. Review "
                                "inputs and try again.").format(args.input_bam))

def _check_input_bam_valid(args, log=None):
    try:
        bamfile = samtools.alignment_file(args.input_bam, 'rb')
        bamfile.close()
    except ValueError:
        raise utils.UsageError(("Specified input [{}] not a valid BAM. Review "
                                "inputs and try again.").format(args.input_bam))

def _check_input_bam_indexed(args, log=None):
    bamfile = samtools.alignment_file(args.input_bam, 'rb')
    try:
        bamfile.fetch()
    except ValueError:
        raise utils.UsageError(("Specified input [{}] is not indexed. Review "
                                "inputs and try again.").format(args.input_bam))
    finally:
        bamfile.close()

def _check_input_bam_not_connor_generated(args, log=None):
    bamfile = samtools.alignment_file(args.input_bam, 'rb')
    header = bamfile.header
    bamfile.close()
    for pg_item in header.get('PG', []):
        if pg_item.get('PN', None) == samtools.CONNOR_PG_PN:
            if args.force:
                msg = ('Specified input [{}] has already been processed with '
                       'Connor (**forcing**)').format(args.input_bam)
                log.warning(msg)
                return
            else:
                msg = ('Specified input [{}] has already been processed with '
                   'Connor. Are you sure? '
                   '(--force to proceed)').format(args.input_bam)
                raise utils.UsageError(msg)

def _check_input_bam_not_empty(args, log=None):
    bamfile = samtools.alignment_file(args.input_bam, 'rb')
    try:
        next(bamfile.fetch())
    except StopIteration:
        raise utils.UsageError(("Specified input [{}] is empty").format(args.input_bam))
    finally:
        bamfile.close()

def _check_input_bam_paired(args, log=None):
    bamfile = samtools.alignment_file(args.input_bam, 'rb')
    try:
        for alignment in bamfile.fetch():
            if alignment.is_paired:
                return
        raise utils.UsageError(("Specified input [{}] does not contain paired reads").format(args.input_bam))
    finally:
        bamfile.close()



def _check_input_bam_barcoded(args, log=None):
    bamfile = samtools.alignment_file(args.input_bam, 'rb')
    softclipped_l = 0
    softclipped_r = 0
    threshold = 0.80
    total_l = 0
    total_r = 0
    sample_size = 1000
    for alignment in itertools.islice(bamfile.fetch(), sample_size):
        if alignment.is_reverse:
            total_r += 1
            if alignment.cigartuples[-1][0] == 4:
                softclipped_r += 1
        else:
            total_l += 1
            if alignment.cigartuples[0][0] == 4:
                softclipped_l += 1
    left_freq = softclipped_l/total_l
    right_freq = softclipped_r/total_r
    if left_freq < threshold or right_freq < threshold:
        raise utils.UsageError(("Specified input [{}] reads do not appear to have barcodes").format(args.input_bam))


_VALIDATIONS = [_check_input_bam_exists,
                _check_input_bam_valid,
                _check_input_bam_indexed,
                _check_input_bam_not_connor_generated,
                _check_input_bam_not_empty,
                _check_input_bam_paired]


def preflight(args, log):
    for validate in _VALIDATIONS:
        validate(args, log)
        