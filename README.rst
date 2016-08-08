======
Connor
======

Command-line tool to deduplicate bam files based on custom, inline barcoding.

.. image:: https://travis-ci.org/umich-brcf-bioinf/Connor.svg?branch=develop
    :target: https://travis-ci.com/umich-brcf-bioinf/Connor
    :alt: Build Status

.. image:: https://codeclimate.com/repos/5793a84516ba097bda009574/badges/28ae96f1f3179a08413e/coverage.svg
   :target: https://codeclimate.com/repos/5793a84516ba097bda009574/coverage
   :alt: Test Coverage

.. image:: https://codeclimate.com/repos/5793a84516ba097bda009574/badges/28ae96f1f3179a08413e/gpa.svg
   :target: https://codeclimate.com/repos/5793a84516ba097bda009574/feed
   :alt: Code Climate


The official repository is at:
https://github.com/umich-brcf-bioinf/Connor

--------
Overview
--------

Connor de-duplicates a barcoded BAM and emits a BAM of consensus alignment pairs
that represent the sequences of original biological molecules. In terms of data
workflow, Connor is similar to Picard MarkDuplicates (ref)).

Sequencing [FASTQ 1/2] -> Aligner [BAM] -> Connor [BAM] -> Variant Detection [VCF]

Consensus output BAM tags:

* X0: filter (rationale explaining why the align was excluded)
* X1: unique identifier for this alignment family
* X2: L~R UMI barcodes for this alignment family; because of fuzzy matching the
  family UMI may be distinct from the UMI of the original alignment
* X3: family size (number of align pairs in this family)
* X4: presence of this tag signals that this alignment would be the template
  for the consensus alignment

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
     -d UMI_DISTANCE_THRESHOLD, --umi_distance_threshold UMI_DISTANCE_THRESHOLD
                           =1 (>=0); UMIs equal to or closer than this Hamming distance will be
                            combined into a single family. Lower threshold make more families with more
                            consistent UMIs; 0 implies UMI must match
                            exactly.

====

Email bfx-connor@umich.edu for support and questions.

UM BRCF Bioinformatics Core
