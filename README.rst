======
Connor
======

Command-line tool to deduplicate reads in bam files based on custom inline barcoding.

.. image:: https://travis-ci.org/umich-brcf-bioinf/Connor.svg?branch=develop
    :target: https://travis-ci.org/umich-brcf-bioinf/Connor
    :alt: Build Status

The official repository is at:
https://github.com/umich-brcf-bioinf/Connor

--------
Overview
--------

TBD

-----------
Quick Start
-----------

TBD

-----------
Connor help
-----------
::
  $ connor --help
  usage: connor input_bam output_bam
  
  Deduplicates BAM file based on custom inline DNA barcodes.
  Emits a new BAM file reduced to a single consensus read for each family of
  original reads.
  
  positional arguments:
    input_bam      path to input BAM
    output_bam     path to output BAM
  
  optional arguments:
    -h, --help     show this help message and exit
    -V, --version  show program's version number and exit

====

Email bfx-connor@umich.edu for support and questions.

UM BRCF Bioinformatics Core
