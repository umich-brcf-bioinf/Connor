import unittest
from connor import connor

class MicroMock(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):
        return 1


class MockLogger(object):
    def __init__(self):
        self._log_calls = []

    def _log(self, msg_format, *args):
        self._log_calls.append((msg_format, args))


class BaseConnorTestCase(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.original_logger = connor._log
        self.mock_logger = MockLogger()
        connor._log = self.mock_logger._log

    def tearDown(self):
        connor._log = self.original_logger
        unittest.TestCase.tearDown(self)
