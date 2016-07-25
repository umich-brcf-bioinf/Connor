from __future__ import print_function, absolute_import, division

from datetime import datetime
import errno
import getpass
import itertools
import logging
import os
import socket
import sys


def zrange(*args):
    try:
        return xrange(*args)
    except NameError:
        return range(*args)

def iter_map(*args):
    try:
        return itertools.imap(*args)
    except AttributeError:
        return map(*args)


class UsageError(Exception):
    """Raised for malformed command or invalid arguments."""
    def __init__(self, msg, *args):
        super(UsageError, self).__init__(msg, *args)



def _makepath(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def _validate_output_file(log_file):
    try:
        _makepath(os.path.dirname(log_file))
        log = open(log_file, "w")
        log.close()
    except Exception:
        raise UsageError(("Connor cannot create specified output file "
                                "[{}]. Review inputs and try again."), log_file)

class Logger(object):
    _DATE_FORMAT = '%Y-%m-%d %H:%M:%S'
    _FILE_LOG_FORMAT = ('%(asctime)s|%(levelname)s|%(start_time)s|%(host)s|'
                        '%(user)s|%(message)s')
    _CONSOLE_LOG_FORMAT = '%(asctime)s|%(levelname)s|%(message)s'

    def __init__(self, args):
        self._verbose = args.verbose
        self._log_filename = args.log_file
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
        self.warning_occurred = False

    def _print(self, level, message, args):
        now = datetime.now().strftime(Logger._DATE_FORMAT)
        print(Logger._CONSOLE_LOG_FORMAT % {'asctime': now,
                                            'levelname': level,
                                            'message': self._format(message,
                                                                    args)},
              file=sys.stderr)
        sys.stderr.flush()

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
        logging.debug(self._format(message, args), extra=self._logging_dict)

    def error(self, message, *args):
        self._print("ERROR", message, args)
        logging.error(self._format(message, args), extra=self._logging_dict)

    def info(self, message, *args):
        self._print("INFO", message, args)
        logging.info(self._format(message, args), extra=self._logging_dict)

    def warning(self, message, *args):
        self._print("WARNING", message, args)
        logging.warning(self._format(message, args), extra=self._logging_dict)
        self.warning_occurred = True
