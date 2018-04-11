======
Connor
======

A command-line tool to deduplicate bam files based on custom, inline barcoding.

.. image:: https://travis-ci.org/umich-brcf-bioinf/Connor.svg?branch=develop
    :target: https://travis-ci.org/umich-brcf-bioinf/Connor
    :alt: Build Status

.. image:: https://codeclimate.com/repos/5793a84516ba097bda009574/badges/28ae96f1f3179a08413e/coverage.svg
    :target: https://codeclimate.com/repos/5793a84516ba097bda009574/coverage
    :alt: Test Coverage

.. image:: https://codeclimate.com/repos/5793a84516ba097bda009574/badges/28ae96f1f3179a08413e/gpa.svg
    :target: https://codeclimate.com/repos/5793a84516ba097bda009574/feed
    :alt: Code Climate

.. image:: https://img.shields.io/badge/license-Apache-blue.svg
    :target: https://pypi.python.org/pypi/connor/
    :alt: License

.. image:: http://img.shields.io/pypi/v/connor.svg
    :target: https://pypi.python.org/pypi/connor/
    :alt: Latest PyPI version

The official repository is at:
https://github.com/umich-brcf-bioinf/Connor

--------
Overview
--------

When analyzing deep-sequence NGS data it is sometimes difficult to distinguish
sequencing and PCR errors from rare variants; as a result some variants may
be missed and some will be identified with an inaccurate variant frequency. To
address this, researchers can attach random barcode sequences during sample
preparation. Upon sequencing, the barcodes act as a signature to trace the set
of PCR amplified molecules back to the original biological molecules of
interest thereby differentiating rare variants in the original molecule from
errors introduced downstream.

Connor accepts a barcoded, paired alignment file (BAM), groups those input
alignments into families, combines each family into a consensus alignment, and
emits the set of deduplicated, consensus alignments (BAM).

   *Connor workflow:*

   *Sequencing [FASTQ 1/2] -> Aligner [BAM] -> Connor [BAM] -> Variant Detection [VCF]*

Connor first groups original alignments into **alignment families** based on their
alignment position and **Universal Molecular Tag (UMT)** barcode. (Connor assumes
the incoming aligned sequences begin with the UMT barcode.) Each family of
alignments is then combined into a single **consensus alignment**; discrepancies
in base-calls and qualities are resolved by majority vote across family members.
By default, smaller families (<3 align pairs) are excluded.

For more information see:

* `QUICKSTART`_ : get started deduplicating barcoded BAMs.

* `INSTALL`_ : alternative ways to install.

* `METHODS`_ : details on UMT barcode structure, suggestions on
  alignment parameters, details on family grouping, and examples.


-----------
Connor help
-----------

::

  $ connor --help
   usage: connor input_bam output_bam

   positional arguments:
     input_bam             path to input BAM
     output_bam            path to deduplicated output BAM

   optional arguments:
     -h, --help            show this help message and exit
     -V, --version         show program's version number and exit
     -v, --verbose         print all log messages to console
     --log_file LOG_FILE   ={output_filename}.log. Path to verbose log file
     --annotated_output_bam ANNOTATED_OUTPUT_BAM
                           path to output BAM containing all original aligns annotated with BAM tags
     -f CONSENSUS_FREQ_THRESHOLD, --consensus_freq_threshold CONSENSUS_FREQ_THRESHOLD
                           =0.6 (0..1.0): Ambiguous base calls at a specific position in a family are
                            transformed to either majority base call, or N if the majority percentage
                            is below this threshold. (Higher threshold results in more Ns in
                            consensus.)
     -s MIN_FAMILY_SIZE_THRESHOLD, --min_family_size_threshold MIN_FAMILY_SIZE_THRESHOLD
                           =3 (>=0): families with count of original reads < threshold are excluded
                            from the deduplicated output. (Higher threshold is more
                            stringent.)
     -d UMT_DISTANCE_THRESHOLD, --umt_distance_threshold UMT_DISTANCE_THRESHOLD
                           =1 (>=0); UMTs equal to or closer than this Hamming distance will be
                            combined into a single family. Lower threshold make more families with more
                            consistent UMTs; 0 implies UMT must match
                            exactly.
     --filter_order {count,name}
                          =count; determines how filters will be ordered in the log
                          results
     --umt_length UMT_LENGTH
                          =6 (>=1); length of UMT

====

Email bfx-connor@umich.edu for support and questions.

UM BRCF Bioinformatics Core

.. _INSTALL: doc/INSTALL.rst
.. _METHODS: doc/METHODS.rst
.. _QUICKSTART : doc/QUICKSTART.rst
