======
Connor
======

A command-line tool to deduplicate bam files based on custom, inline barcoding.

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

See `METHODS`_ for more details on UMT barcode structure, suggestions on
alignment parameters, details on family grouping, and examples.

-----------
Quick Start
-----------

1. **Install Connor**

The following command will install connor and its dependencies. The command may
take several moments to start showing progress and several minutes to complete.
See `INSTALL`_ for more info on requirements and more install options.
::
  $ pip install git+https://github.com/umich-brcf-bioinf/Connor
  ...
  Successfully installed Connor-0.3
  
  $ connor --help
  usage: connor input_bam output_bam
  positional arguments:
  ...
  v0.3
 
2. **Download the example data:**
::
  $ git clone https://github.com/umich-brcf-bioinf/Connor sample_data
  ...

3. **Run Connor:**

This command will read PIK3CA-original.bam and produce PIK3CA-deduped.bam (with
BAM index) and a log file in your working directory.

::
  $ connor sample_data/examples/PIK3CA-original.bam PIK3CA-deduped.bam
  2016-08-11 17:08:07|INFO|connor begins
  2016-08-11 17:08:07|INFO|logging to [PIK3CA-deduped.bam.log]
  2016-08-11 17:08:07|INFO|reading input bam [sample_data/examples/PIK3CA-original.bam]
  ...
  2016-08-11 17:09:06|INFO|1016/3482 (29.18%) families were excluded because the original read count < 3
  2016-08-11 17:09:06|INFO|61937 original pairs were deduplicated to 2466 families (overall dedup rate 96.02%)
  2016-08-11 17:09:06|INFO|2466 families written to [PIK3CA-deduped.bam]
  2016-08-11 17:09:06|INFO|connor complete

  $ ls -1 PIK3CA*
  PIK3CA-deduped.bam
  PIK3CA-deduped.bam.bai
  PIK3CA-deduped.bam.log

From the log above, you can see that the original alignments resulted in 2466
deduplicated families (pairs).

4. **Reviewing output file**

If you have samtools installed, you can examine the difference between original
and deduplicated bams:

::
  $ samtools flagstat Connor/examples/PIK3CA-original.bam
  158401 + 0 in total (QC-passed reads + QC-failed reads)
  722 + 0 secondary
  0 + 0 supplementary
  0 + 0 duplicates
  158401 + 0 mapped (100.00%:nan%)
  157679 + 0 paired in sequencing
  78817 + 0 read1
  78862 + 0 read2
  155079 + 0 properly paired (98.35%:nan%)
  157362 + 0 with itself and mate mapped
  317 + 0 singletons (0.20%:nan%)
  2145 + 0 with mate mapped to a different chr
  1991 + 0 with mate mapped to a different chr (mapQ>=5)
  
  $ samtools flagstat PIK3CA-deduped.bam
  4932 + 0 in total (QC-passed reads + QC-failed reads)
  0 + 0 secondary
  0 + 0 supplementary
  0 + 0 duplicates
  4932 + 0 mapped (100.00%:nan%)
  4932 + 0 paired in sequencing
  2466 + 0 read1
  2466 + 0 read2
  4932 + 0 properly paired (100.00%:nan%)
  4932 + 0 with itself and mate mapped
  0 + 0 singletons (0.00%:nan%)
  0 + 0 with mate mapped to a different chr
  0 + 0 with mate mapped to a different chr (mapQ>=5)

Note that 158401 original alignments were deduplicated to 4932 (2466 pairs).

5. **Reviewing a consensus alignment**

Each consensus alignment has a set of custom tags that provide details
on the family of original alignments.
::
  $ samtools view PIK3CA-deduped.bam | head -1 | tr '\t' '\n'
  NS500501:108:HMKNKBGXX:1:13205:18985:5894
  163
  chr3
  ...
  X1:i:175
  X2:Z:ATGGAT~AAGACC
  X3:i:41

The documentation for these tags is in the SAM/BAM header and excerpted here:

* X1: unique identifier (integer) for this alignment family
* X2: Left~Right UMT barcodes for this alignment family; because of fuzzy matching the
  family UMT may be distinct from the UMT of the original alignment
* X3: family size (number of align pairs in this family)

Interpreting the tag definitions with the alignment above, the consensus
alignment **175** (X1) represents **41** original alignment pairs (X3) whose
alignment position matched exactly and left-right UMT barcodes matched
**ATGGAT-AAGACC** (X2).

-----------
Connor help
-----------

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
     -d UMI_DISTANCE_THRESHOLD, --umi_distance_threshold UMI_DISTANCE_THRESHOLD
                           =1 (>=0); UMIs equal to or closer than this Hamming distance will be
                            combined into a single family. Lower threshold make more families with more
                            consistent UMIs; 0 implies UMI must match
                            exactly.

====

Email bfx-connor@umich.edu for support and questions.

UM BRCF Bioinformatics Core

.. _METHODS: METHODS.rst
.. _INSTALL: INSTALL.rst


