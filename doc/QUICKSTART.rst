-----------
Quick Start
-----------

1. **Install Python**

Launch a terminal window, and check that Python (2.7 or higher) and pip are installed.
::
  $ python --version
  $ pip --version

If Python or pip is not available, see instructions for installing here:
 * https://www.python.org/downloads/
 * https://packaging.python.org/installing/#install-pip-setuptools-and-wheel

2. **Install Connor**

This command will install connor and its dependencies; it may take a moment to start
showing progress and should take several minutes to complete. See `INSTALL`_ for more
info on requirements and more install options.
::
  $ pip install connor --user
  ...
  Successfully installed Connor-0.4

Run this line (and/or add it to your ~/.bashrc file) and confirm connor is in your path:
::
  $ PATH=~/.local/bin:$PATH
  $ which connor
  <some_dir>/connor
  $ connor --help
  usage: connor input_bam output_bam
  positional arguments:
  ...
  v0.5
 
3. **Download the example data:**
::
  $ git clone https://github.com/umich-brcf-bioinf/Connor sample_data
  ...

4. **Run Connor:**

This command will read PIK3CA-original.bam and produce PIK3CA-deduped.bam (with
BAM index) and a log file in your working directory.
::
  $ connor sample_data/examples/PIK3CA-original.bam PIK3CA-deduped.bam
  2016-08-11 17:08:07|INFO|connor begins
  2016-08-11 17:08:07|INFO|logging to [PIK3CA-deduped.bam.log]
  2016-08-11 17:08:07|INFO|reading input bam [sample_data/examples/PIK3CA-original.bam]
  ...
  2016-08-11 17:08:29|INFO|8.69% (17937/206509) alignments unplaced or discarded
  2016-09-11 15:08:29|INFO|families discarded: 16.03% (1518/9470) family too small (<3)
  2016-09-11 15:08:29|INFO|91.31% (188572/206509) alignments included in 7952 families
  2016-09-11 15:08:29|INFO|95.78% deduplication rate (1 - 7952 families/188572 included alignments)
  2016-09-11 15:08:30|INFO|sorting/indexing [examples/PIK3CA-deduped.bam]
  2016-08-11 15:08:30|INFO|connor complete  (23 seconds, 146mb peak memory)

  $ ls -1 PIK3CA*
  PIK3CA-deduped.bam
  PIK3CA-deduped.bam.bai
  PIK3CA-deduped.bam.log

From the log above, you can see that the original alignments resulted in 7952
deduplicated families (pairs).

5. **Reviewing output file**

You don't need samtools installed to use Connor, but if you have samtools installed,
you can examine the difference between original and deduplicated bams:
::
  $ samtools flagstat sample_data/examples/PIK3CA-original.bam
  206509 + 0 in total (QC-passed reads + QC-failed reads)
  59 + 0 secondary
  0 + 0 supplementary
  0 + 0 duplicates
  205705 + 0 mapped (99.61%:nan%)
  206450 + 0 paired in sequencing
  103175 + 0 read1
  103275 + 0 read2
  203974 + 0 properly paired (98.80%:nan%)
  204842 + 0 with itself and mate mapped
  804 + 0 singletons (0.39%:nan%)
  677 + 0 with mate mapped to a different chr
  664 + 0 with mate mapped to a different chr (mapQ>=5)
  
  $ samtools flagstat PIK3CA-deduped.bam
  15904 + 0 in total (QC-passed reads + QC-failed reads)
  0 + 0 secondary
  0 + 0 supplementary
  0 + 0 duplicates
  15904 + 0 mapped (100.00%:nan%)
  15904 + 0 paired in sequencing
  7952 + 0 read1
  7952 + 0 read2
  15904 + 0 properly paired (100.00%:nan%)
  15904 + 0 with itself and mate mapped
  0 + 0 singletons (0.00%:nan%)
  0 + 0 with mate mapped to a different chr
  0 + 0 with mate mapped to a different chr (mapQ>=5)

Note that 206509 original alignments were deduplicated to 15904 (7952 pairs).

6. **Reviewing a consensus alignment**

Connor adds a set of custom tags to each consensus alignment that provide details
on the family of original alignment. Here is an excerpt of the first alignment:
::
  $ samtools view PIK3CA-deduped.bam | head -1 | tr '\t' '\n'
  HWI-D00143:749:HM5YFBCXX:2:1112:3541:48875
  99
  chr3
  178873584
  ... lines omitted ...
  X1:Z:178873584~178873660
  X2:Z:14S76M10S~8S76M16S
  X3:i:0
  X4:Z:GAAAGT~CTTCGT
  X5:i:5


Note the BAM tags listed on the last few lines above (X1-X5). The documentation for these
tags is in the SAM/BAM header and summarized here:

* X1: leftmost~rightmost matched pair positions
* X2: L~R CIGARs
* X3: unique identifier for this alignment family
* X4: L~R UMT barcodes for this alignment family; because of fuzzy matching the 
  family UMT may be distinct from the UMT of the original alignment
* X5: family size (number of align pairs in this family)

Interpreting the tag definitions with the alignment above, this consensus
alignment represents **5** original alignment pairs (from tag X5 above) whose
alignment leftmost and rightmost positions matched **178873584~178873660**
(from tag X1) and left-right UMT barcodes matched **ATGGAT-AAGACC** (from
tag X4).

.. _METHODS: METHODS.rst
.. _INSTALL: INSTALL.rst
