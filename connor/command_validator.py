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
    #pylint: disable=redefined-builtin
    import itertools.izip as iter_zip
except ImportError:
    iter_zip = zip     #pylint: disable=invalid-name
import os
import connor.utils as utils
import connor.samtools as samtools

_SAMPLE_SIZE = 1000

def _balanced_strand_gen(aligns, total_aligns):
    '''given collection of random {aligns}, returns alternating pos/neg
    strands up to {total_aligns}; will return smaller of pos/neg collection
    if less than total_aligns'''
    predicate=lambda align: align.is_reverse
    gen1, gen2 = itertools.tee((predicate(align), align) for align in aligns)
    pos_strand_gen = (align for pred, align in gen1 if not pred)
    neg_strand_gen = (align for pred, align in gen2 if pred)
    combined_gen = iter_zip(pos_strand_gen, neg_strand_gen)
    interleaved_gen = (align for aligns in combined_gen for align in aligns)
    return itertools.islice(interleaved_gen, total_aligns)

def _log_force_or_raise(args, log, msg):
    if args.force:
        log.warning(msg + ' (**forcing**)')
    else:
        raise utils.UsageError(msg + ' Are you sure? (--force to proceed)')

def _sample_bamfile(input_bam, extractor_function):
    stats = {'forward': [], 'reverse': []}
    bamfile = samtools.alignment_file(input_bam, 'rb')
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
           'reads.').format(args.input_bam)
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

    softclipped = _sample_bamfile(args.input_bam, is_edge_softclipped)
    percent_true = lambda stats, strand : sum(stats[strand])/len(stats[strand])
    if _below_threshold(softclipped, percent_true, SOFTCLIP_THRESHOLD):
        msg = ('Specified input [{}] reads do not appear to have '
                   'barcodes.').format(args.input_bam)
        _log_force_or_raise(args, log, msg)

def _check_input_bam_consistent_length(args, log=None):
    #pylint: disable=invalid-name
    CONSISTENT_LENGTH_THRESHOLD = 0.90

    def freq_of_most_common_length(seq_lengths, strand):
        lengths = seq_lengths[strand]
        counter = Counter(lengths)
        most_common_length_count = counter.most_common(1)[0]
        count = most_common_length_count[1]
        return count / len(lengths)

    seq_lengths = _sample_bamfile(args.input_bam,
                                 lambda a: len(a.query_sequence))
    if _below_threshold(seq_lengths,
                        freq_of_most_common_length,
                        CONSISTENT_LENGTH_THRESHOLD):
        msg = ('Specified input [{}] reads appear to have inconsistent '
               'sequence lengths.').format(args.input_bam)
        _log_force_or_raise(args, log, msg)

_VALIDATIONS = [_check_input_bam_exists,
                _check_input_bam_valid,
                _check_input_bam_indexed,
                _check_input_bam_not_deduped,
                _check_input_bam_not_empty,
                _check_input_bam_paired,
                _check_input_bam_consistent_length]


def preflight(args, log):
    for validate in _VALIDATIONS:
        validate(args, log)
        