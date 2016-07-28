#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
import unittest
from collections import defaultdict
from collections import OrderedDict
import connor.utils as utils

class MicroMock(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):
        return 1


class MockLogger(object):
    def __init__(self):
        self._log_calls = defaultdict(list)

    def error(self, message, *args):
        self._log_calls["ERROR"].append(self._format(message, args))

    def warning(self, message, *args):
        self._log_calls["WARNING"].append(self._format(message, args))

    def info(self, message, *args):
        self._log_calls["INFO"].append(self._format(message, args))

    def debug(self, message, *args):
        self._log_calls["DEBUG"].append(self._format(message, args))

    @staticmethod
    def _format(message, args):
        return message.format(*[i for i in args])


class BaseConnorTestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.mock_logger = MockLogger()

    def tearDown(self):
        unittest.TestCase.tearDown(self)


class FilteredGeneratorTest(BaseConnorTestCase):
    def test_filter_passesAllThroughWhenNoFilters(self):
        base_collection = [1,2,3,4,5]
        filters = {}

        generator = utils.FilteredGenerator(filters)
        actual_collection  = [x for x in generator.filter(base_collection)]

        self.assertEqual(base_collection, actual_collection)
        self.assertEqual(5, generator.total_included)
        self.assertEqual(0, generator.total_excluded)
        self.assertEqual(0, len(generator.filter_stats))

    def test_filter_singleFilter(self):
        filters = {'div by 2' : lambda x: x % 2 == 0}
        generator = utils.FilteredGenerator(filters)

        base_collection = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        actual_collection  = [x for x in generator.filter(base_collection)]

        self.assertEqual([1,3,5,7,9], actual_collection)
        self.assertEqual(5, generator.total_included)
        self.assertEqual(5, generator.total_excluded)
        self.assertEqual(5, generator.filter_stats['div by 2'])
        self.assertEqual(1, len(generator.filter_stats))

    def test_filter_multipleFilters(self):
        filters = {'div by 2': lambda x: x % 2 == 0,
                   'div by 5': lambda x: x % 5 == 0}
        generator = utils.FilteredGenerator(filters)

        base_collection = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        actual_collection  = [x for x in generator.filter(base_collection)]

        self.assertEqual([1,3,7,9], actual_collection)
        self.assertEqual(4, generator.total_included)
        self.assertEqual(6, generator.total_excluded)
        self.assertEqual(4, generator.filter_stats['div by 2'])
        self.assertEqual(1, generator.filter_stats['div by 5'])
        self.assertEqual(1, generator.filter_stats['div by 2;div by 5'])
        self.assertEqual(3, len(generator.filter_stats))

    def test_filter_multipleFilterAreSortedByName(self):
        filters = {'div by 5': lambda x: x % 5 == 0,
                   'div by 2': lambda x: x % 2 == 0}
        generator = utils.FilteredGenerator(filters)

        base_collection = [1, 10]
        actual_collection  = [x for x in generator.filter(base_collection)]

        self.assertEqual([1], actual_collection)
        self.assertEqual(1, generator.filter_stats['div by 2;div by 5'])

    def test_filter_stats_immutable(self):
        generator = utils.FilteredGenerator({})
        try:
            generator.filter_stats = {'foo': 42}
            self.assertEqual(True, False, 'filter_stats should be immutable')
        except AttributeError:
            self.assertEqual(True, True, 'filter_stats is immutable')

        generator.filter_stats['foo'] = 42
        base_collection = [1, 2 ,3]
        actual_collection  = [x for x in generator.filter(base_collection)]

        self.assertEqual(base_collection, actual_collection)
        self.assertEqual(0, len(generator.filter_stats))

    def test_filter_stats_ordersResultsByCountBreakingTiesByName(self):
        generator = utils.FilteredGenerator({})
        generator._filter_stats = {'C': 1,
                                   'B;C':2,
                                   'B': 2,
                                   'E': 3,
                                   'A': 3,
                                   'D': 4}
        actual_filter_stats = generator.filter_stats
        self.assertEqual(6, len(actual_filter_stats))
        expected_filter_stats = OrderedDict([('D', 4),
                                             ('A', 3),
                                             ('E' ,3),
                                             ('B', 2),
                                             ('B;C', 2),
                                             ('C', 1)])
        self.assertEqual(expected_filter_stats,
                         actual_filter_stats)


    def test_filter_returnsEmptyIfBaseEmpty(self):
        generator = utils.FilteredGenerator({})
        base_collection = []
        actual_collection  = [x for x in generator.filter(base_collection)]

        self.assertEqual([], actual_collection)
        self.assertEqual(0, len(generator.filter_stats))
