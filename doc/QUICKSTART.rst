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
  v0.4
 
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

5. **Reviewing output file**

You don't need samtools installed to use Connor, but if you have samtools installed,
you can examine the difference between original and deduplicated bams:
::
  $ samtools flagstat sample_data/examples/PIK3CA-original.bam
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

6. **Reviewing a consensus alignment**

Connor adds a set of custom tags to each consensus alignment that provide details
on the family of original alignment. Here is an excerpt of the first alignment:
::
  $ samtools view PIK3CA-deduped.bam | head -1 | tr '\t' '\n'
  NS500501:108:HMKNKBGXX:1:13205:18985:5894
  163
  chr3
  . . .
  X1:i:175
  X2:Z:ATGGAT~AAGACC
  X3:i:41

Note the BAM tags listed on the last few lines above (X1-X3). The documentation for these
tags is in the SAM/BAM header and summarized here:

* X1: unique identifier (integer) for this alignment family
* X2: Left~Right UMT barcodes for this alignment family; because of fuzzy matching the
  family UMT may be distinct from the UMT of the original alignment
* X3: family size (number of align pairs in this family)

Interpreting the tag definitions with the alignment above, this consensus
alignment represents **41** original alignment pairs (from tag X3 above) whose
alignment position matched exactly and left-right UMT barcodes matched
**ATGGAT-AAGACC** (from tag X2).

.. _METHODS: METHODS.rst
.. _INSTALL: INSTALL.rst
