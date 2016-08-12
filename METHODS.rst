Background
==========

Researchers are increasingly searching for lower allele frequency
variants in diverse genomic populations using deep coverage
next-generation sequencing data(refs). Errors in (a) PCR duplication and
(b) instrument base calling confound attempts to reduce false positive
variants and assign true biological alternate allele frequency for rare
variants. In response, a number of researchers have introduced custom
DNA barcoding (refs). By adding a molecule-specific barcode sequence,
these methods enable downstream data transformations that can
meaningfully aggregate sequences of PCR duplicates into sequence of the
original biological molecules and disambiguate PCR and sequencing errors
from low-frequency variants.

Connor de-duplicates a barcoded BAM and emits a BAM of consensus
alignment pairs that represent the sequences of original biological
molecules. In terms of data workflow, Connor is similar to
position-based deduplication (e.g. Picard MarkDuplicates (ref)).
However, deduplication based solely on position disproportionally
discards original molecules that fall in small target regions with
ultra-deep coverage (where positional overlap is more frequent) and
typically chooses a consensus sequences that match the reference over
sequences containing variants. In contrast, Connor combines sequences
where the alignment structure and molecular barcodes match, creating a
consensus sequences that model the population of original molecules.

Workflow
========

    Sequencing [FASTQ 1/2] -> Aligner [BAM] -> Connor [BAM] -> Variant
    Detection [VCF]

Connor assumes the input BAM files are barcoded using the Rubicon
ThruPLEX R044 Preparation Kit (ref). In particular Connor assumes each
query sequence begins with a Universal Molecular Tag (UMT) in-line
barcode sequence (6 nt), followed by a non-random stem sequence (8-11
nt), followed by the target sequence region. The aligner should preserve
the leading UMT barcode and stem sequences and mark those areas as "soft
clipped" in the BAM CIGAR field to indicate they did not match the
reference; this is default behavior for some aligners (e.g.
Burrows-Wheeler Aligner, BWA).

+--------------------+-----------------------+-------------------+---------------------------+
|     **region**     |     **UMT barcode**   |     **stem**      |     **target sequence**   |
+====================+=======================+===================+===========================+
|     **sequence**   |     ACTGTT            |     GTAGCTCA      |     GTTGAGACACAT...       |
+--------------------+-----------------------+-------------------+---------------------------+
|     **CIGAR**      |     soft clipped      |     matched ref   |
+--------------------+-----------------------+-------------------+---------------------------+

Because a correct UMT and consistent alignment structure are integral to
Connor’s ability to accurately deduplicate, we encourage analysts to
avoid any manipulations between sequencing and alignment and
deduplication. Given a set of PCR replicates from a single original
molecule, examples of problematic manipulations include:

-  end trimming can shave base calls off the front of the sequence (the
   location of the UMT barcode) which would prevent affected reads from
   being correctly grouped

-  pre-alignment quality trimming may shave reads off the end of the
   sequence; this would create distinct CIGAR values obscuring the
   original alignment structure

Connor has been tested with the following aligners (using default
parameters except where noted):

-  BWA v. 0.7.12

-  Bowtie2 v. 2.2.4 (-local mode)

-  DNASTAR SeqMan NGen v.13.0.2.2 (disable clipping and deduplication)

-  Hisat2 v.2.0.4

-  Novoalign v.3.04.06

Overview of UMT barcode deduplication method
============================================

1. Discard alignments which could not be mapped, are not properly
   paired, have low mapping quality (<1), or are missing CIGAR value.

2. Discard alignments whose pair partner is missing

3. Group together alignment pairs which share the same stem-template
   edge coordinates into “position families”

4. Based on left and right UMT, subdivide each “position family” into
   “UMT families”.

   a. Extract left + right (combined 12-mer) UMT barcodes and sort by
      frequency into candidate families

   b. For each alignment pair, loop over candidate family (in descending
      popularity) comparing the alignment UMT with the candidate family

      -  Exact match across 12-mers is considered a match

      -  Exact match of left *or* right UMT (6-mer) is considered a
         match

      -  *Inexact* match (within a user-defined Hamming distance,
         default=1) of left or right UMT is considered a match

   c. Each alignment pair will either match an existing family, or found
      its own family

5. Within each UMT family, establish the majority CIGAR across
   alignments; discard alignment pairs with non-majority CIGARs

6. Discard UMT families with fewer than a user-defined number of
   original pairs (default=3 pairs)

7. For each UMT family, collapse the set of original alignment pairs
   into a single consensus alignment pair.

-  Consensus sequence is determined by majority vote at each position in
   the base call sequence

-  Any position with less than a user-defined percent majority
   (default=60%) results in an N at that position

-  Consensus quality is determined by majority vote at each position in
   the quality sequence

-  Consensus CIGAR is the majority cigar for that UMT family

Figures
=======

**Figure 1: Alignment pairs are grouped into position families**

|fig1_position_families| 

Reads are grouped into position families based on their
alignment to the reference. Each bar represents a paired alignment
positioned on the reference. The stem-target boundaries of the target
sequences define the groups (beginning of left target, end of the right
target) [a,b] Unmapped, unpaired, or low-quality reads are removed.

**Figure 2: Position families are subdivided into UMT families**

|fig2_umt_families|

Position families are further divided into UMT families by UMT barcode.
Each bar represents a paired alignment; distinct UMT barcodes are
numbered and shaded. The collection of the left is the alignments from
the position family; on the right those same alignments are partitioned
into three UMT families. Rare anomalies in tagging and PCR processes can
erroneously split families. To prevent spurious family splits, an
alignment will be admitted to a UMT family if either left or right UMT
barcode is similar (within specified hamming distance). Several match
types are noted: [a] are exact left-right matches, [b] is exact one-side
match, [c] is inexact one-side match. Small families (<3 alignment
pairs) are discarded [d].

**
**

**Figure 3: UMT family of alignments are combined into a consensus**

|fig3_consensus|

The set of alignment pairs in a UMT family (top) are combined into a
single consensus pair (bottom) by majority vote. Each bar represents a
paired alignment; distinct UMT barcodes are numbered and shaded;
mismatches (candidate variants) in the target sequence region are
highlighted. Within a UMT family, only structurally identical alignments
(i.e. matching CIGAR values) can be combined; alignment pairs with
minority CIGAR values are discarded [a]. Consensus alignment preserves
the majority UMT left and right barcodes and stem. In the target
sequence region, the base calls are tallied for each position and the
majority base call becomes the consensus base call [b,c]. If the
majority base call is less than the consensus sequence threshold (60% by
default), the base call is replaced by “N” indicating ambiguity [d].

References
==========

TBD

.. |fig1_position_families| image:: doc/images/fig1_position_families.png
.. |fig2_umt_families| image:: doc/images/fig2_umt_families.png
.. |fig3_consensus| image:: doc/images/fig3_consensus.png
