from __future__ import print_function, absolute_import, division
import math

def build_family_handlers(args,
                          consensus_writer,
                          annotated_writer,
                          logger):
    handlers = [_FamilySizeStatHandler(logger),
                _MatchStatHandler(args, logger),
                _WriteAnnotatedAlignsHandler(annotated_writer),
                _WriteConsensusHandler(consensus_writer)]
    return handlers

class _WriteConsensusHandler(object):
    def __init__(self, consensus_writer):
        self._writer = consensus_writer
    def handle(self, family):
        if not family.filter_value:
            self._writer.write(family,
                               family.consensus,
                               family.consensus.left)
            self._writer.write(family,
                               family.consensus,
                               family.consensus.right)

    def end(self):
        pass


class _WriteAnnotatedAlignsHandler(object):
    def __init__(self, writer):
        self._writer = writer

    def handle(self, family):
        for align_pair in sorted(family.align_pairs,
                                 key=lambda a: a.query_name):
            self._writer.write(family,
                               align_pair,
                               align_pair.left)
            self._writer.write(family,
                               align_pair,
                               align_pair.right)

    def end(self):
        pass


class _FamilySizeStatHandler(object):
    def __init__(self, logger):
        self.collection = []
        self.min = None
        self.quartile_1 = None
        self.median = None
        self.mean = None
        self.quartile_3 = None
        self.max = None
        self._log = logger

    def handle(self, tag_family):
        self.collection.append(len(tag_family.align_pairs))

    @staticmethod
    def _percentile(collection, percent):
        if not collection:
            return None
        fractional_index = (len(collection)-1) * percent
        floor_index = int(math.floor(fractional_index))
        ceiling_index = int(math.ceil(fractional_index))
        index = int(fractional_index)
        if floor_index == ceiling_index:
            value = collection[index]
        else:
            lower = collection[floor_index] * (ceiling_index - fractional_index)
            upper = collection[ceiling_index] * (fractional_index - floor_index)
            value = lower + upper
        return value

    @property
    def summary(self):
        return (self.min,
                self.quartile_1,
                self.median,
                self.mean,
                self.quartile_3,
                self.max)

    #TODO: cgates: add guard for empty collection
    def end(self):
        percentile = _FamilySizeStatHandler._percentile
        self.collection.sort()
        self.min = self.collection[0]
        self.quartile_1 = percentile(self.collection, 0.25)
        self.median = percentile(self.collection, 0.50)
        self.quartile_3 = percentile(self.collection, 0.75)
        self.max = self.collection[-1]
        self.mean = sum(self.collection) / len(self.collection)

        self._log.debug(('family_stat|family size distribution (original pair '
                   'counts: min, 1Q, median, mean, 3Q, max): {}'),
                  ', '.join(map(lambda x: str(round(x, 2)), self.summary)))


class _MatchStatHandler(object):
    def __init__(self, args, logger):
        self._log = logger
        self.hamming_threshold = args.umt_distance_threshold
        self.total_inexact_match_count = 0
        self.total_pair_count = 0

    def handle(self, tag_family):
        self.total_inexact_match_count += tag_family.inexact_match_count
        self.total_pair_count += len(tag_family.align_pairs)

    def end(self):
        exact_count = self.total_pair_count - self.total_inexact_match_count
        self._log.debug(('{:.2f}% ({}/{}) original pairs matched '
                         'UMT exactly'),
                        100 * (1 - self.percent_inexact_match),
                        exact_count,
                        self.total_pair_count)

        self._log.debug(('{:.2f}% ({}/{}) original pairs matched '
                      'by Hamming distance threshold (<={}) on '
                      'left or right UMT '),
                     100 * self.percent_inexact_match,
                     self.total_inexact_match_count,
                     self.total_pair_count,
                     self.hamming_threshold)

    @property
    def percent_inexact_match(self):
        return self.total_inexact_match_count/self.total_pair_count
