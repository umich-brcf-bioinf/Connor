from __future__ import print_function, absolute_import, division
from collections import defaultdict
import pandas as pd
import connor.samtools as samtools

def build_family_handlers(args,
                          consensus_writer,
                          annotated_writer,
                          logger):
    handlers = [_FamilySizeStatHandler(logger),
                _MatchStatHandler(args, logger),
                _CigarMinorityStatHandler(logger),
                _CigarStatHandler(logger),
                _WriteConsensusHandler(consensus_writer)]
    if annotated_writer != samtools.AlignWriter.NULL:
        handlers.append(_WriteAnnotatedAlignsHandler(annotated_writer))
    return handlers

class _WriteConsensusHandler(object):
    def __init__(self, consensus_writer):
        self._writer = consensus_writer
    def handle(self, family):
        if not family.filter_value:
            self._writer.write(family,
                               family.consensus.left)
            self._writer.write(family,
                               family.consensus.right)

    def end(self):
        pass

# class _WriteFamilyHandler(object):
#     def __init__(self, args, consensus_writer, logger):
#         self._writer = consensus_writer
#         self._min_family_size_threshold = args.min_family_size_threshold
#         self._log = logger
#         self.included_family_count = 0
#         self.excluded_family_count = 0
#         self.total_alignment_count = 0
# 
#     def handle(self, tag_family):
#         self.total_alignment_count += len(tag_family.align_pairs)
#         #TODO: cgates: This check should be done using the filter field
#         if len(tag_family.align_pairs) >= self._min_family_size_threshold:
#             consensus = tag_family.consensus
#             self._writer.write(tag_family,
#                                consensus.left)
#             self._writer.write(tag_family,
#                                consensus.right)
#             self.included_family_count += 1
#         else:
#             self.excluded_family_count += 1
# 
#     def end(self):
#         total_family_count = self.excluded_family_count + \
#                              self.included_family_count
#         self._log.info(('{}/{} ({:.2f}%) families were excluded because the '
#                    'original read count < {}'),
#                   self.excluded_family_count,
#                   total_family_count,
#              100 * self.excluded_family_count / total_family_count,
#              self._min_family_size_threshold)
#         dedup_percent = 100 * \
#                 (1 - (self.included_family_count / self.total_alignment_count))
#         self._log.info(('{} original pairs were deduplicated to {} families '
#                    '(overall dedup rate {:.2f}%)'),
#                   self.total_alignment_count,
#                   self.included_family_count,
#                   dedup_percent)
#         self._log.info('{} families written to [{}]',
#                   self.included_family_count,
#                   self._writer.bam_file_path)

class _WriteAnnotatedAlignsHandler(object):
    def __init__(self, writer):
        self._writer = writer

    def handle(self, family):
        for align_pair in family.align_pairs:
            self._writer.write(family, align_pair.left)
            self._writer.write(family, align_pair.right)

    def end(self):
        pass


class _BaseTukeyStatHandler(object):
    def __init__(self):
        self.collection = []
        self.min = None
        self.quartile_1 = None
        self.median = None
        self.quartile_3 = None
        self.max = None

    def get_family_statistic(self, tag_family):
        pass

    def handle(self, tag_family):
        stat = self.get_family_statistic(tag_family)
        if stat is not None:
            self.collection.append(stat)

    @property
    def summary(self):
        return (self.min,
                self.quartile_1,
                self.median,
                self.quartile_3,
                self.max)

    def end(self):
        summary = pd.Series(self.collection).describe()
        self.min = summary['min']
        self.max = summary['max']
        self.median = summary['50%']
        self.quartile_1 = summary['25%']
        self.quartile_3 = summary['75%']


class _FamilySizeStatHandler(_BaseTukeyStatHandler):
    def __init__(self, logger):
        super(_FamilySizeStatHandler, self).__init__()
        self._log = logger

    def get_family_statistic(self, tag_family):
        return len(tag_family.align_pairs)

    def end(self):
        super(_FamilySizeStatHandler, self).end()
        self._log.debug(('family_stat|family size distribution (original pair '
                   'counts: min, 1Q, median, 3Q, max): {}'),
                  ', '.join(map(str, self.summary)))


class _CigarMinorityStatHandler(_BaseTukeyStatHandler):
    def __init__(self, logger):
        super(_CigarMinorityStatHandler, self).__init__()
        self._log = logger

    def get_family_statistic(self, tag_family):
        stat = tag_family.minority_cigar_percentage
        if stat:
            return stat
        else:
            return None

    def end(self):
        super(_CigarMinorityStatHandler, self).end()
        self._log.debug(('family_stat|cigar|family distribution of minority '
                   'CIGAR percentages (min, 1Q, median, 3Q, max): {}'),
                  ', '.join(map(lambda x: str(round(x,2)), self.summary)))

#TODO: (cgates): Split into tukey and other stats
class _CigarStatHandler(object):
    def __init__(self, logger):
        self._log = logger
        self.distinct_cigar_counts = []
        self.min = None
        self.quartile_1 = None
        self.median = None
        self.quartile_3 = None
        self.max = None
        self.total_input_alignment_count = 0
        self.total_alignment_count = 0
        self.total_family_count = 0
        self.families_cigar_counts = defaultdict(int)
        self.family_cigar_minority_percentage_counts = defaultdict(list)


    def handle(self, tag_family):
        self.distinct_cigar_counts.append(tag_family.distinct_cigar_count)
        self.total_input_alignment_count += len(tag_family.align_pairs)
        self.total_alignment_count += tag_family.included_pair_count
        self.total_family_count += 1
        self.families_cigar_counts[tag_family.distinct_cigar_count] += 1
        self.family_cigar_minority_percentage_counts[tag_family.distinct_cigar_count].append(tag_family.minority_cigar_percentage)


    @property
    def percent_deduplication(self):
        return 1 - (self.total_family_count / self.total_input_alignment_count)


    @property
    def total_excluded_alignments(self):
        return self.total_input_alignment_count - self.total_alignment_count

    @property
    def percent_excluded_alignments(self):
        return self.total_excluded_alignments / self.total_input_alignment_count

    @property
    def summary(self):
        return (self.min,
                self.quartile_1,
                self.median,
                self.quartile_3,
                self.max)

    def end(self):
        summary = pd.Series(self.distinct_cigar_counts).describe()
        self.min = summary['min']
        self.max = summary['max']
        self.median = summary['50%']
        self.quartile_1 = summary['25%']
        self.quartile_3 = summary['75%']

        self._log.debug(('family_stat|cigar|family distribution of distinct'
                         ' CIGAR counts (min, 1Q, median, 3Q, max): {}'),
                        ', '.join(map(str, self.summary)))
        ordered_cigar_counts = sorted(self.families_cigar_counts.items(),
                                      key = lambda x: -1 * int(x[1]))
        for num_cigars, freq in ordered_cigar_counts:
            summary = pd.Series(self.family_cigar_minority_percentage_counts[num_cigars]).describe()
            self._log.debug(('family_stat|cigar|{}/{} ({:.2f}%) families had '
                             '{} CIGAR: minor % distrib '
                             '{:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}'),
                            freq,
                            self.total_family_count,
                            100 * freq / self.total_family_count,
                            num_cigars,
                            summary['min'],
                            summary['25%'],
                            summary['50%'],
                            summary['75%'],
                            summary['max'])
        self._log.debug(('family_stat|cigar|{}/{} ({:.2f}%) pairs were '
                         'excluded as minority CIGAR'),
                        self.total_excluded_alignments,
                        self.total_input_alignment_count,
                        100 * self.percent_excluded_alignments)
        self._log.debug(('family_stat|{} original pairs (of majority CIGAR) '
                         'were deduplicated to {} families '
                         '(majority CIGAR dedup rate {:.2f}%)'),
                        self.total_alignment_count,
                        self.total_family_count,
                        100 * self.percent_deduplication)


class _MatchStatHandler(object):
    def __init__(self, args, logger):
        self._log = logger
        self.hamming_threshold = args.umi_distance_threshold
        self.total_inexact_match_count = 0
        self.total_pair_count = 0

    def handle(self, tag_family):
        self.total_inexact_match_count += tag_family.inexact_match_count
        self.total_pair_count += len(tag_family.align_pairs)

    def end(self):
        exact_count = self.total_pair_count - self.total_inexact_match_count
        self._log.debug(('family_stat|{}/{} ({:.2f}%) original pairs matched '
                         'UMI exactly'),
                        exact_count,
                        self.total_pair_count,
                        100 * (1 - self.percent_inexact_match))

        self._log.debug(('family_stat|{}/{} ({:.2f}%) original pairs matched '
                      'by Hamming distance threshold (<={}) on '
                      'left or right UMI '),
                     self.total_inexact_match_count,
                     self.total_pair_count,
                     100 * self.percent_inexact_match,
                     self.hamming_threshold)

    @property
    def percent_inexact_match(self):
        return self.total_inexact_match_count/self.total_pair_count
