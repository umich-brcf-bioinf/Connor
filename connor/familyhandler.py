from __future__ import print_function, absolute_import, division
from collections import defaultdict
import pandas as pd
import connor.samtools as samtools

#TODO: (cgates): switch all handlers to family-at-at-time handling
def build_family_handlers(args, logger):
    handlers = [_FamilySizeStatHandler(logger),
            _MatchStatHandler(args, logger),
            _CigarMinorityStatHandler(logger),
            _CigarStatHandler(logger),
            _WriteFamilyHandler(args, logger)]
    if args.output_excluded_alignments:
        handlers.append(_WriteExcludedReadsHandler(args, logger))
    return handlers

class _WriteFamilyHandler(object):
    def __init__(self, args, logger):
        self._output_filename = args.output_bam
        input_bamfile = samtools.alignment_file(args.input_bam, 'rb')
        self._output_bamfile = samtools.alignment_file(self._output_filename,
                                                       'wb',
                                                       template=input_bamfile)
        input_bamfile.close()
        self._min_family_size_threshold = args.min_family_size_threshold
        self._log = logger
        self.included_family_count = 0
        self.excluded_family_count = 0
        self.total_alignment_count = 0

    def handle(self, tag_families):
        for tag_family in tag_families:
            self.total_alignment_count += len(tag_family.alignments)
            if len(tag_family.alignments) >= self._min_family_size_threshold:
                consensus_pair = tag_family.consensus
                self._output_bamfile.write(consensus_pair.left_alignment)
                self._output_bamfile.write(consensus_pair.right_alignment)
                self.included_family_count += 1
            else:
                self.excluded_family_count += 1

    def end(self):
        total_family_count = self.excluded_family_count + \
                             self.included_family_count
        self._log.info(('{}/{} ({:.2f}%) families were excluded because the '
                   'original read count < {}'),
                  self.excluded_family_count,
                  total_family_count,
             100 * self.excluded_family_count / total_family_count,
             self._min_family_size_threshold)
        dedup_percent = 100 * (1 - (self.included_family_count / self.total_alignment_count))
        self._log.info(('{} original pairs were deduplicated to {} families '
                   '(overall dedup rate {:.2f}%)'),
                  self.total_alignment_count,
                  self.included_family_count,
                  dedup_percent)
        self._log.info('{} families written to [{}]',
                  self.included_family_count,
                  self._output_filename)
        self._output_bamfile.close()
        self._log.info('sorting and indexing [{}]', self._output_filename)
        samtools.sort_and_index_bam(self._output_filename)


class _WriteExcludedReadsHandler(object):
    def __init__(self, args, logger):
        self._output_filename = args.output_excluded_alignments
        input_bamfile = samtools.alignment_file(args.input_bam, 'rb')
        self._output_bamfile = samtools.alignment_file(self._output_filename,
                                                       'wb',
                                                       template=input_bamfile)
        input_bamfile.close()
        self._min_family_size_threshold = args.min_family_size_threshold
        self._log = logger
        self.total_alignment_count = 0

    def handle(self, tag_families):
        for tag_family in tag_families:
            self.total_alignment_count += len(tag_family.excluded_alignments)
            self.total_alignment_count += 1
            self.tag_reads(tag_family, tag_family.consensus, included=1)
            self._output_bamfile.write(tag_family.consensus.left_alignment)
            self._output_bamfile.write(tag_family.consensus.right_alignment)
            for pair in tag_family.excluded_alignments:
                self.tag_reads(tag_family, pair, included=0)
                self._output_bamfile.write(pair.left_alignment)
                self._output_bamfile.write(pair.right_alignment)

    def tag_reads(self, tag_family, original_pair, included):
        x0 = tag_family.umi_sequence
        x1 = "{0}|{1}".format(tag_family.umi[0], tag_family.umi[1])
        x2 = "{0},{1}".format(original_pair.left_alignment.reference_start + 1,
                              original_pair.right_alignment.reference_end)
        x3 = len(tag_family.alignments)
        x4 = str(x3 >= self._min_family_size_threshold)
        original_pair.left_alignment.set_tag("X0", x0, "i")
        original_pair.left_alignment.set_tag("X1", x1, "Z")
        original_pair.left_alignment.set_tag("X2", x2, "Z")
        original_pair.left_alignment.set_tag("X3", x3, "i")
        original_pair.left_alignment.set_tag("X4", x4, "Z")
        original_pair.left_alignment.set_tag("X5", included, "i")

        original_pair.right_alignment.set_tag("X0", x0, "i")
        original_pair.right_alignment.set_tag("X1", x1, "Z")
        original_pair.right_alignment.set_tag("X2", x2, "Z")
        original_pair.right_alignment.set_tag("X3", x3, "i")
        original_pair.right_alignment.set_tag("X4", x4, "Z")
        original_pair.right_alignment.set_tag("X5", included, "i")

    def end(self):
        self._log.info('{} tagged original pairs written to [{}]',
                  self.total_alignment_count,
                  self._output_filename)
        self._output_bamfile.close()
        self._log.info('sorting and indexing [{}]', self._output_filename)
        samtools.sort_and_index_bam(self._output_filename)

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

    def handle(self, tag_families):
        for tag_family in tag_families:
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
        return len(tag_family.alignments)

    def end(self):
        super(_FamilySizeStatHandler, self).end()
        self._log.debug(('family stat|family size distribution (original pair '
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
        self._log.debug(('family stat|cigar|family distribution of minority '
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


    def handle(self, tag_families):
        for tag_family in tag_families:
            self.distinct_cigar_counts.append(tag_family.distinct_cigar_count)
            self.total_input_alignment_count += tag_family.input_alignment_count
            self.total_alignment_count += len(tag_family.alignments)
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

        self._log.debug(('family stat|cigar|family distribution of distinct'
                         ' CIGAR counts (min, 1Q, median, 3Q, max): {}'),
                        ', '.join(map(str, self.summary)))
        ordered_cigar_counts = sorted(self.families_cigar_counts.items(),
                                      key = lambda x: -1 * int(x[1]))
        for num_cigars, freq in ordered_cigar_counts:
            summary = pd.Series(self.family_cigar_minority_percentage_counts[num_cigars]).describe()
            self._log.debug(('family stat|cigar|{}/{} ({:.2f}%) families had '
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
        self._log.debug(('family stat|cigar|{}/{} ({:.2f}%) pairs were '
                         'excluded as minority CIGAR'),
                        self.total_excluded_alignments,
                        self.total_input_alignment_count,
                        100 * self.percent_excluded_alignments)
        self._log.debug(('family stat|{} original pairs (of majority CIGAR) '
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

    def handle(self, tag_families):
        for tag_family in tag_families:
            self.total_inexact_match_count += tag_family.inexact_match_count
            self.total_pair_count += len(tag_family.alignments)

    def end(self):
        exact_match_count = self.total_pair_count - self.total_inexact_match_count
        self._log.debug(('family stat|{}/{} ({:.2f}%) original pairs matched '
                         'UMI exactly'),
                        exact_match_count,
                        self.total_pair_count,
                        100 * (1 - self.percent_inexact_match))

        self._log.debug(('family stat|{}/{} ({:.2f}%) original pairs matched '
                      'by Hamming distance threshold (<={}) on '
                      'left or right UMI '),
                     self.total_inexact_match_count,
                     self.total_pair_count,
                     100 * self.percent_inexact_match,
                     self.hamming_threshold)

    @property
    def percent_inexact_match(self):
        return self.total_inexact_match_count/self.total_pair_count
