======
Connor
======

Command-line tool to deduplicate reads in bam files based on custom inline barcoding.

.. image:: https://travis-ci.org/umich-brcf-bioinf/Connor.svg?branch=develop
    :target: https://travis-ci.com/umich-brcf-bioinf/Connor
    :alt: Build Status

.. image:: https://coveralls.io/repos/github/umich-brcf-bioinf/Connor/badge.svg?branch=develop
    :target: https://coveralls.io/github/umich-brcf-bioinf/Connor?branch=develop
    :alt: Coverage Status


The official repository is at:
https://github.com/umich-brcf-bioinf/Connor

--------
Overview
--------

TBD

Consensus output BAM tags:
  X0: UMI sequence: a numeric identifier assigned to a consensus sequence
  X1: UMI tag: left|right inline barcodes
  X2: edge positions: start of the left alignment, end of the right alignment
  X3: count of original read pairs that were deduplicated into the consensus
  
-----------
Quick Start
-----------

1. **Install Connor (see INSTALL.rst):**
::
  $ pip install git+https://github.com/umich-brcf-bioinf/Connor

2. **Get the examples directory:**
::
  $ git clone https://github.com/umich-brcf-bioinf/Connor

3. **Run Connor:**
::
  $ connor Connor/examples/PIK3CA-original.bam PIK3CA-deduped.bam

This will read PIK3CA-original.bam and produce PIK3CA-deduped.bam (in your
working directory). 


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
    -f CONSENSUS_FREQ_THRESHOLD, --consensus_freq_threshold CONSENSUS_FREQ_THRESHOLD
                        =0.6 (0..1.0): Ambiguous base calls at a specific position in a family are 
                        transformed to either majority base call, or N if the majority percentage 
                        is below this threshold. (Higher threshold results in more Ns in consensus.)
    -s MIN_FAMILY_SIZE_THRESHOLD, --min_family_size_threshold MIN_FAMILY_SIZE_THRESHOLD
                        =3 (>=0): families with count of original reads < threshold are excluded
                        from the deduplicated output. (Higher threshold is more stringent.)
    -d UMI_DISTANCE_THRESHOLD, --umi_distance_threshold UMI_DISTANCE_THRESHOLD
                        =1 (>=0); UMIs equal to or closer than this Hamming distance will be 
                        combined into a single family. Lower threshold make more families with more 
                        consistent UMIs; 0 implies UMI must match exactly.

====

Email bfx-connor@umich.edu for support and questions.

UM BRCF Bioinformatics Core
