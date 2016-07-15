#pylint: disable=line-too-long
import os
from setuptools import find_packages, setup
import connor

def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as filename:
        return filename.read()

setup(name='Connor',
      version = connor.__version__,
      description=('Command-line tool to deduplicate reads in bam files based '
                   'on custom inline barcoding.'),
      long_description=(read('README.rst') + '\n\n' +
                        read('CHANGELOG.rst') + '\n\n' +
                        read('AUTHORS.rst')),
      url='https://github.com/umich-brcf-bioinf/Connor',
      author='University of Michigan Bioinformatics Core',
      author_email='bfx-connor@umich.edu',
      license='Apache',
      packages=find_packages(exclude=['test*']),
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: Apache Software License',
                   'Operating System :: Unix',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'],
      keywords='bioinformatic exome-seq DNA-seq BAM',
      setup_requires=['cython'],
      install_requires=['pysam>=0.8.4,<0.9.1', 'pandas>=0.18'],
      entry_points={'console_scripts': ['connor=connor.connor:main']},
      test_suite='nose.collector',
      tests_require=['nose', 'pysam', 'testfixtures'],
      zip_safe=False)
