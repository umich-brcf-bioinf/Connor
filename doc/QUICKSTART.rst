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

.. _METHODS: METHODS.rst
.. _INSTALL: INSTALL.rst
