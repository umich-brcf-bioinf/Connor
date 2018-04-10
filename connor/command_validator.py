'''Validates command preconditions.
Specifically checks that the command, arguments, and environment
(e.g. input/output directories or files) are consistent and plausible.
Each validation function evaluates a specific precondition raising a
UsageException if things look problematic or logging a warning if the condition
can be forced.
'''
from __future__ import print_function, absolute_import, division
from collections import Counter
import itertools
try:
    from itertools import izip as iter_zip
except ImportError:
    #pylint: disable=invalid-name
    iter_zip = zip
import os
import connor.consam.writers as writers
import connor.consam.pysamwrapper as pysamwrapper
import connor.utils as utils

_SAMPLE_SIZE = 1000

def _balanced_strand_gen(base_aligns, limit):
    '''given collection of random {aligns}, returns alternating forward/reverse
    strands up to {total_aligns}; will return smaller of forward/reverse
    collection if less than total_aligns'''
    predicate=lambda align: align.is_reverse
    gen1, gen2 = itertools.tee((predicate(align), align) for align in base_aligns)
    for_strand_gen = (align for pred, align in gen1 if not pred)
    rev_strand_gen = (align for pred, align in gen2 if pred)
    combined_gen = iter_zip(for_strand_gen, rev_strand_gen)
    interleaved_gen = (align for for_rev_aligns in combined_gen for align in for_rev_aligns)
    return itertools.islice(interleaved_gen, limit)

def _check_stat_below_threshold(args,
                                extractor_function,
                                freq_function,
                                threshold,
                                msg,
                                log):
    stats = _sample_bamfile(args.input_bam,
                            extractor_function)
    if _below_threshold(stats,
                    freq_function,
                    threshold):
        _log_force_or_raise(args, log, msg)

def _log_force_or_raise(args, log, msg):
    if args.force:
        log.warning(msg + ' (**forcing**)')
    else:
        raise utils.UsageError(msg + ' Are you sure? (--force to proceed)')

def _sample_bamfile(input_bam, extractor_function):
    stats = {'forward': [], 'reverse': []}
    bamfile = pysamwrapper.alignment_file(input_bam, 'rb')
    try:
        for align in _balanced_strand_gen(bamfile.fetch(), _SAMPLE_SIZE):
            stats[_strand(align)].append(extractor_function(align))
    finally:
        bamfile.close()
    return stats

def _below_threshold(strand_stats,
                     freq_function,
                     threshold):
    forward_freq = freq_function(strand_stats, 'forward')
    reverse_freq = freq_function(strand_stats, 'reverse')
    return min(forward_freq, reverse_freq) < threshold

def _check_input_bam_exists(args, log=None): #pylint: disable=unused-argument
    if not os.path.exists(args.input_bam):
        raise utils.UsageError(("Specified input [{}] does not exist. Review "
                                "inputs and try again.").format(args.input_bam))

def _check_input_bam_valid(args, log=None): #pylint: disable=unused-argument
    try:
        bamfile = pysamwrapper.alignment_file(args.input_bam, 'rb')
        bamfile.close()
    except ValueError:
        raise utils.UsageError(("Specified input [{}] not a valid BAM. Review "
                                "inputs and try again.").format(args.input_bam))

def _check_input_bam_indexed(args, log=None): #pylint: disable=unused-argument
    bamfile = pysamwrapper.alignment_file(args.input_bam, 'rb')
    try:
        bamfile.fetch()
    except ValueError:
        raise utils.UsageError(("Specified input [{}] is not indexed. Review "
                                "inputs and try again.").format(args.input_bam))
    finally:
        bamfile.close()

def _check_input_bam_not_deduped(args, log=None):
    bamfile = pysamwrapper.alignment_file(args.input_bam, 'rb')
    header_dict = pysamwrapper.get_header_dict(bamfile)
    bamfile.close()
    names = set([pg_item.get('PN', None) for pg_item in header_dict.get('PG', [])])
    if writers.CONNOR_PG_PN in names:
        msg = ('Specified input [{}] has already been processed with '
               'Connor.').format(args.input_bam)
        _log_force_or_raise(args, log, msg)

def _check_input_bam_not_empty(args, log=None): #pylint: disable=unused-argument
    bamfile = pysamwrapper.alignment_file(args.input_bam, 'rb')
    try:
        next(bamfile.fetch())
    except StopIteration:
        msg = "Specified input [{}] is empty"
        raise utils.UsageError(msg.format(args.input_bam))
    finally:
        bamfile.close()

def _check_input_bam_paired(args, log=None):
    bamfile = pysamwrapper.alignment_file(args.input_bam, 'rb')
    try:
        for alignment in itertools.islice(bamfile.fetch(), _SAMPLE_SIZE):
            if alignment.is_paired:
                return
    finally:
        bamfile.close()
    msg = ('Specified input [{}] does not appear to contain paired '
           'reads.').format(args.input_bam)
    _log_force_or_raise(args, log, msg)

def _check_input_bam_properly_paired(args, log=None):
    bamfile = pysamwrapper.alignment_file(args.input_bam, 'rb')
    try:
        for alignment in itertools.islice(bamfile.fetch(), _SAMPLE_SIZE):
            if alignment.is_proper_pair:
                return
    finally:
        bamfile.close()
    msg = ('Specified input [{}] does not appear to contain any properly paired '
           'alignments.').format(args.input_bam)
    _log_force_or_raise(args, log, msg)

def _strand(align):
    return 'reverse' if align.is_reverse else 'forward'


def _check_input_bam_barcoded(args, log=None):
    #pylint: disable=invalid-name
    SOFTCLIP_OP = 4
    SOFTCLIP_THRESHOLD = 0.80

    def is_edge_softclipped(align):
        edge_index = -1 if align.is_reverse else 0
        return int(align.cigartuples[edge_index][0] == SOFTCLIP_OP)

    percent_true = lambda stats, strand : sum(stats[strand])/len(stats[strand])
    msg = ('Specified input [{}] reads do not appear to have '
               'barcodes.').format(args.input_bam)
    _check_stat_below_threshold(args,
                                is_edge_softclipped,
                                percent_true,
                                SOFTCLIP_THRESHOLD,
                                msg,
                                log)

def _check_input_bam_consistent_length(args, log=None):
    #pylint: disable=invalid-name
    CONSISTENT_LENGTH_THRESHOLD = 0.90

    def freq_of_most_common_length(seq_lengths, strand):
        lengths = seq_lengths[strand]
        counter = Counter(lengths)
        most_common_length_count = counter.most_common(1)[0]
        count = most_common_length_count[1]
        return count / len(lengths)

    msg = ('Specified input [{}] reads appear to have inconsistent '
           'sequence lengths.').format(args.input_bam)
    _check_stat_below_threshold(args,
                                lambda a: len(a.query_sequence),
                                freq_of_most_common_length,
                                CONSISTENT_LENGTH_THRESHOLD,
                                msg,
                                log)

def _check_overwrite_output(args, log=None):
    output_collisions = []
    if os.path.exists(args.output_bam):
        output_collisions.append(args.output_bam)
    if args.annotated_output_bam and os.path.exists(args.annotated_output_bam):
        output_collisions.append(args.annotated_output_bam)
    if output_collisions:
        output_collision_str = ', '.join(output_collisions)
        msg = 'One or more outputs [{}] exist.'.format(output_collision_str)
        _log_force_or_raise(args, log, msg)

_VALIDATIONS = [_check_input_bam_exists,
                _check_input_bam_valid,
                _check_input_bam_indexed,
                _check_input_bam_not_deduped,
                _check_input_bam_not_empty,
                _check_input_bam_paired,
                _check_input_bam_properly_paired,
                _check_input_bam_consistent_length,
                _check_overwrite_output]


def preflight(args, log):
    for validate in _VALIDATIONS:
        validate(args, log)
