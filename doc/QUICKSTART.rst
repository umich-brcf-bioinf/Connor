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
  2016-08-25 11:22:38|INFO|connor begins (v0.3)
  2016-08-25 11:22:38|INFO|logging to [PIK3CA-deduped.bam.log]
  2016-08-25 11:22:38|INFO|reading input bam [sample_data/examples/PIK3CA-original.bam]
  ... (some lines omitted) ...
  2016-08-25 11:23:04|INFO|16.00% (1515/9468) families discarded: family too small (<3)
  2016-08-25 11:23:04|INFO|91.32% (188578/206509) alignments included in 7953 families
  2016-08-25 11:23:04|INFO|95.78% deduplication rate (1 - 7953 families/188578 included alignments)
  2016-08-25 11:23:09|INFO|sorting/indexing [PIK3CA-deduped.bam]
  2016-08-25 11:23:10|INFO|connor complete (32 seconds, 194mb peak memory).

  $ ls -1 PIK3CA*
  PIK3CA-deduped.bam
  PIK3CA-deduped.bam.bai
  PIK3CA-deduped.bam.log

From the log above, you can see that the original alignments resulted in 7953
deduplicated families (pairs).

4. **Reviewing output file**

If you have samtools installed, you can examine the difference between original
and deduplicated bams:
::
  $ samtools flagstat Connor/examples/PIK3CA-original.bam
  206509 + 0 in total (QC-passed reads + QC-failed reads)
  59 + 0 secondary
  0 + 0 supplementary
  0 + 0 duplicates
  205705 + 0 mapped (99.61%:-nan%)
  206450 + 0 paired in sequencing
  103175 + 0 read1
  103275 + 0 read2
  203974 + 0 properly paired (98.80%:-nan%)
  204842 + 0 with itself and mate mapped
  804 + 0 singletons (0.39%:-nan%)
  677 + 0 with mate mapped to a different chr
  664 + 0 with mate mapped to a different chr (mapQ>=5)
  
  $ samtools flagstat PIK3CA-deduped.bam
  15906 + 0 in total (QC-passed reads + QC-failed reads)
  0 + 0 secondary
  0 + 0 supplementary
  0 + 0 duplicates
  15906 + 0 mapped (100.00%:-nan%)
  15906 + 0 paired in sequencing
  7953 + 0 read1
  7953 + 0 read2
  15906 + 0 properly paired (100.00%:-nan%)
  15906 + 0 with itself and mate mapped
  0 + 0 singletons (0.00%:-nan%)
  0 + 0 with mate mapped to a different chr
  0 + 0 with mate mapped to a different chr (mapQ>=5)

Note that 206509 original alignments were deduplicated to 15906 (7953 pairs).

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

  HWI-D00143:749:HM5YFBCXX:2:1112:3541:48875
  99
  chr3
  ... (some lines omitted) ...
  X1:i:0
  X2:Z:GAAAGT~CTTCGT
  X3:i:5
  
The documentation for these tags is in the SAM/BAM header and excerpted here:

* X1: unique identifier (integer) for this alignment family
* X2: Left~Right UMT barcodes for this alignment family; because of fuzzy matching the
  family UMT may be distinct from the UMT of the original alignment
* X3: family size (number of align pairs in this family)

Interpreting the tag definitions with the alignment above, this consensus
alignment represents **5** original alignment pairs (X3) whose
alignment position matched exactly and left-right UMT barcodes matched
**ATGGAT-AAGACC** (X2).

.. _METHODS: METHODS.rst
.. _INSTALL: INSTALL.rst
