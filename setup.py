from __future__ import print_function
import os
import platform
import sys

from setuptools import find_packages, setup

import connor

_REQUIRED_PYTHON_VERSION = (2, 7)

def check_python_version():
    if sys.version_info < _REQUIRED_PYTHON_VERSION:
        msg_format = '''
Problem: Python v{0}.{1} or above is required but you are using v{2}.
Please install a supported version of Python and try again.\
'''
        message = msg_format.format(_REQUIRED_PYTHON_VERSION[0],
                                    _REQUIRED_PYTHON_VERSION[1],
                                    platform.python_version())
        print(message, file=sys.stderr)
        sys.exit(1)

def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as filename:
        return filename.read()

check_python_version()

setup(name='Connor',
      version = connor.__version__,
      description=('Command-line tool to deduplicate reads in bam files based '
                   'on custom inline barcoding.'),
      long_description=(read('README.rst') + '\n\n' +
                        read('doc/CHANGELOG.rst') + '\n\n' +
                        read('doc/AUTHORS.rst')),
      url='https://github.com/umich-brcf-bioinf/Connor',
      author='University of Michigan Bioinformatics Core',
      author_email='bfx-connor@umich.edu',
      license='Apache',
      packages=find_packages(exclude=['test*']),
      classifiers=['Development Status :: 5 - Production/Stable',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: Apache Software License',
                   'Operating System :: Unix',
                   'Operating System :: MacOS',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'],
      keywords='bioinformatic exome-seq DNA-seq BAM',
      setup_requires=['cython'],
      install_requires=['pysam>=0.8.4,<0.9.1', 'sortedcontainers>=1.5.3'],
      entry_points={'console_scripts': ['connor=connor.connor:main']},
      test_suite='nose.collector',
      tests_require=['nose', 'pysam', 'testfixtures'],
      zip_safe=False)
