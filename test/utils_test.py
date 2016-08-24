#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=unused-argument

from argparse import Namespace
from collections import defaultdict
import unittest
import connor.utils as utils


class MicroMock(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):
        return 1


class MockLogger(object):
    def __init__(self, *args):
        self.warning_occurred = False
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

    def ok(self):
        #pylint: disable=redundant-unittest-assert
        self.assertTrue(True)

class FilteredGeneratorTest(BaseConnorTestCase):
    def test_filter_passesAllThroughWhenNoFilters(self):
        base_collection = [1,2,3,4,5]
        filters = {}

        generator = utils.FilteredGenerator(filters)
        actual_collection  = [x for x in generator.filter(base_collection)]

        self.assertEqual([(1, None),
                          (2, None),
                          (3, None),
                          (4, None),
                          (5, None)], actual_collection)

    def test_filter_singleFilter(self):
        filters = {'div by 2' : lambda x: x % 2 == 0}
        generator = utils.FilteredGenerator(filters)

        base_collection = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        actual_collection  = [x for x in generator.filter(base_collection)]

        self.assertEqual([(1, None),
                          (2, 'div by 2'),
                          (3, None),
                          (4, 'div by 2'),
                          (5, None),
                          (6, 'div by 2'),
                          (7, None),
                          (8, 'div by 2'),
                          (9, None),
                          (10, 'div by 2')], actual_collection)

    def test_filter_multipleFilters(self):
        filters = {'div by 2': lambda x: x % 2 == 0,
                   'div by 5': lambda x: x % 5 == 0}
        generator = utils.FilteredGenerator(filters)

        base_collection = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        actual_collection  = [x for x in generator.filter(base_collection)]

        self.assertEqual([(1, None),
                          (2, 'div by 2'),
                          (3, None),
                          (4, 'div by 2'),
                          (5, 'div by 5'),
                          (6, 'div by 2'),
                          (7, None),
                          (8, 'div by 2'),
                          (9, None),
                          (10, 'div by 2;div by 5')],
                         actual_collection)

    def test_filter_multipleFilterAreSortedByName(self):
        filters = {'div by 5': lambda x: x % 5 == 0,
                   'div by 2': lambda x: x % 2 == 0}
        generator = utils.FilteredGenerator(filters)

        base_collection = [1, 10]
        actual_collection  = [x for x in generator.filter(base_collection)]

        self.assertEqual([(1, None),
                          (10, 'div by 2;div by 5')], actual_collection)

    def test_filter_returnsEmptyIfBaseEmpty(self):
        generator = utils.FilteredGenerator({})
        base_collection = []
        actual_collection  = [x for x in generator.filter(base_collection)]

        self.assertEqual([], actual_collection)

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


class MockBaseLogger(object):
    def __init__(self):
        self.lines = []
    def info(self, line, extra=None): #pylint: disable=unused-argument
        self.lines.append(line)
    def debug(self, line, extra=None):
        self.info(line, extra)
    def error(self, line, extra=None):
        self.info(line, extra)
    def warning(self, line, extra=None):
        self.info(line, extra)


class LoggerTestCase(BaseConnorTestCase):
    def test_consoleFormat(self):
        args = Namespace(log_file="test.log",
                         verbose=None)
        console_stream = StringIO()
        logger = utils.Logger(args, console_stream)
        logger._file_logger = MockBaseLogger()
        logger.info("info [{}]", "i")
        console_lines = console_stream.getvalue().strip().split("\n")
        self.assertEqual(1, len(console_lines))
        fields = console_lines[0].split("|")
        self.assertEqual(3, len(fields))
        self.assertEqual(fields[1], "INFO")
        self.assertEqual(fields[2], "info [i]")


    def test_notVerbose(self):
        args = Namespace(log_file="test.log",
                         verbose=None)
        console_stream = StringIO()
        file_logger = MockBaseLogger()
        logger = utils.Logger(args, console_stream)
        logger._file_logger = file_logger
        logger.debug("debug [{}]", "d")
        logger.info("info [{}]", "i")
        logger.warning("warning [{}]", "w")
        logger.error("error [{}]", "e")

        console_lines = console_stream.getvalue().strip().split("\n")
        self.assertEqual(3, len(console_lines))
        console_lines[0].endswith("info [i]")
        console_lines[1].endswith("warning [w]")
        console_lines[2].endswith("error [e]")
        self.assertEqual(4, len(file_logger.lines))
        file_logger.lines[0].endswith("debug [d]")
        file_logger.lines[1].endswith("info [i]")
        file_logger.lines[2].endswith("warning [w]")
        file_logger.lines[3].endswith("error [e]")

    def test_consoleVerbose(self):
        args = Namespace(log_file="test.log",
                         verbose=True)
        console_stream = StringIO()
        logger = utils.Logger(args, console_stream)
        file_logger = MockBaseLogger()
        logger._file_logger = file_logger
        logger.debug("debug [{}]", "d")
        logger.info("info [{}]", "i")
        logger.warning("warning [{}]", "w")
        logger.error("error [{}]", "e")

        console_lines = console_stream.getvalue().strip().split("\n")
        self.assertEqual(4, len(console_lines))
        console_lines[0].endswith("debug [d]")
        console_lines[1].endswith("info [i]")
        console_lines[2].endswith("warning [w]")
        console_lines[3].endswith("error [e]")
        self.assertEqual(4, len(file_logger.lines))
        file_logger.lines[0].endswith("debug [d]")
        file_logger.lines[1].endswith("info [i]")
        file_logger.lines[2].endswith("warning [w]")
        file_logger.lines[3].endswith("error [e]")


    def test_warning_setsWarningOccurred(self):
        args = Namespace(log_file="test.log",
                         verbose=True)
        console_stream = StringIO()
        logger = utils.Logger(args, console_stream)
        logger._file_logger = MockBaseLogger()

        self.assertEqual(False, logger.warning_occurred)
        logger.warning("oops")
        self.assertEqual(True, logger.warning_occurred)
