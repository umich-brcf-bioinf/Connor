Changelog
=========

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

