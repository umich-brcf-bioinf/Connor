Changelog
=========

0.5.x (MM/DD/YYYY)
------------------
- Extended supported python and pysam versions
- Adjusted to avoid performance problem when processing extremely deep pileups
- Adjusted so that when no families pass filters show warning instead of
  error message (thanks to ccario83 for upvoting this fix)

0.5 (9/13/2016)
---------------
- Filters now exclude supplemental alignments
- Added BAM tags to show pair positions and CIGAR values
- Reduced required memory and improved performance

0.4 (8/26/2016)
---------------
- Added input/command validations
- Added annotated bam option
- Revised QUICKSTART, METHODS
- Added PG line in BAM header
- Improved logging of filtered aligns and progress
- Removed some logged stats to focus logging results
- Removed dependency on pandas/numpy
- Moderate performance (speed) improvements in calculating consensus sequence
- Switched consensus quality to be the max mapping quality

0.3 (8/8/2016)
---------------
- Added filters to exclude low quality, unmapped, or unpaired alignments
- Revised BAM tags; documented BAM tags in BAM header
- Extended logging to write to file and console
- Adjusted to make deterministic in Py3/Py2

0.2 (7/15/2016)
---------------
- Bugfix: connor was mangling left hand side of right hand consensus reads
- Fuzzy grouping of pairs into families based on left or right UMI match
- Fuzzy grouping of pairs into families based on UMI within Hamming distance
- Command line args for hamming distince, consensus threshold, min orig reads
- Extended logging to assist in overall diagnostics
- Generate additional file of alignments excluded from consensus (diagnostic)
- Added UMI sequence tag (X0)

0.1 (6/17/2016)
---------------
- Initial development release
- Partitions raw reads into consensus families
