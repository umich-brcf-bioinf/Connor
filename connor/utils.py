from __future__ import print_function, absolute_import, division

from datetime import datetime
import getpass
try:
    from itertools import map as iter_map
except ImportError:
    iter_map = map
import logging
import socket
import sys


class UsageError(Exception):
    '''Raised for malformed command or invalid arguments.'''
    def __init__(self, msg, *args):
        super(UsageError, self).__init__(msg, *args)

class CountingGenerator(object):
    def __init__(self):
        self.item_count = 0

    def count(self, generator):
        for i in generator:
            self.item_count += 1
            yield i

class FilteredGenerator(object):
    '''Applies filters to a base collection yielding the item and its filter'''
    def __init__(self, filter_dict):
        '''
        Args:
            filter_dict (dict): key = filter name, value = function that
                that accepts an item and returns true is that item should
                be excluded. For example: {"div by 2": lambda x: x % 2 == 0}
        '''
        self._filters = sorted(filter_dict.items(), key=lambda x: x[0])

    def filter(self, base_collection):
        '''Yields subset of base_collection/generator based on filters.'''
        for item in base_collection:
            excluded = []
            for (name, exclude) in self._filters:
                if exclude(item):
                    excluded.append(name)
            if excluded:
                filter_value = ";".join(excluded)
            else:
                filter_value = None
            yield item, filter_value


class Logger(object):
    _DATE_FORMAT = '%Y-%m-%d %H:%M:%S'
    _FILE_LOG_FORMAT = ('%(asctime)s|%(levelname)s|%(start_time)s|%(host)s|'
                        '%(user)s|%(message)s')
    _CONSOLE_LOG_FORMAT = '%(asctime)s|%(levelname)s|%(message)s'

    def __init__(self, args, console_stream=None):
        self._verbose = args.verbose
        self._log_filename = args.log_file
        if console_stream:
            self._console_stream = console_stream
        else:
            self._console_stream = sys.stderr
        user = getpass.getuser()
        host = socket.gethostname()
        start_time = datetime.now().strftime(Logger._DATE_FORMAT)
        self._logging_dict = {'user': user,
                              'host': host,
                              'start_time' : start_time}
        logging.basicConfig(format=Logger._FILE_LOG_FORMAT,
                            level="DEBUG",
                            datefmt=Logger._DATE_FORMAT,
                            filename=self._log_filename)
        self._file_logger = logging
        self.warning_occurred = False

    def _print(self, level, message, args):
        now = datetime.now().strftime(Logger._DATE_FORMAT)
        print(Logger._CONSOLE_LOG_FORMAT % {'asctime': now,
                                            'levelname': level,
                                            'message': self._format(message,
                                                                    args)},
              file=self._console_stream)
        self._console_stream.flush()

    @staticmethod
    def _format(message, args):
        try:
            log_message = message.format(*[i for i in args])
        except IndexError as err:
            log_message = ("Malformed log message ({}: {})"
                           "|{}|{}").format(type(err).__name__,
                                            err,
                                            message,
                                            [str(i) for i in args])
        return log_message

    def debug(self, message, *args):
        if self._verbose:
            self._print("DEBUG", message, args)
        self._file_logger.debug(self._format(message, args),
                               extra=self._logging_dict)

    def _log(self, msg_type, method, message, *args):
        self._print(msg_type, message, args)
        method(self._format(message, args),
               extra=self._logging_dict)

    def error(self, message, *args):
        self._log("ERROR", self._file_logger.error, message, *args)

    def info(self, message, *args):
        self._log("INFO", self._file_logger.info, message, *args)

    def warning(self, message, *args):
        self._log("WARNING", self._file_logger.warning, message, *args)
        self.warning_occurred = True

def sort_dict(key_counts):
    return sorted(key_counts.items(), key=lambda x: (-1 * x[1], x[0]))
