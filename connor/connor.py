'''
Created on Jun 3, 2016

@author: pulintz, cgates
'''
from __future__ import print_function, absolute_import, division
from collections import defaultdict
import pysam

class LightweightAlignment(object):
    """ Minimal info from PySam.AlignedSegment"""
    def __init__(self, aligned_segment):
        self.name = aligned_segment.query_name
        chrom = aligned_segment.reference_name
        pos1 = aligned_segment.reference_start
        pos2 = aligned_segment.next_reference_start
        if pos1 < pos2:
            self.key = (chrom, pos1, pos2)
        else:
            self.key = (chrom, pos2, pos1)

def _build_alignment_family_dict(lw_aligns):
    af_dict = defaultdict(set)
    for lwa in lw_aligns:
        af_dict[lwa.key].add(lwa.name)
    return af_dict

def _build_read_families(aligned_segments,coord_read_name_dict):
    family_dict = defaultdict(set)
    for aseg in aligned_segments:
        key = LightweightAlignment(aseg).key
        family_dict[key].add(aseg)
        if (2*len(coord_read_name_dict[key])) == len(family_dict[key]):
            yield family_dict.pop(key)

def _build_consensus_pair(alignments):
    start_alignment = None
    for alignment in alignments:
        if not start_alignment:
            start_alignment = alignment
        elif alignment.query_name == start_alignment.query_name:
            return (start_alignment, alignment)


def main(input_bam, output_bam):
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    lw_aligns = [LightweightAlignment(align) for align in bamfile.fetch()]
    align_family_dict = _build_alignment_family_dict(lw_aligns)
    bamfile.close()
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    outfile = pysam.AlignmentFile(output_bam, "wb", template=bamfile)
    for family in _build_read_families(bamfile.fetch(), align_family_dict):
        read1, read2 = _build_consensus_pair(family)
        outfile.write(read1)
        outfile.write(read2)
    outfile.close()
    bamfile.close()

if __name__ == '__main__':
    input_bam="/Volumes/MyPassport/Data/Rubicon/CU1/BAM/EGFR-ENSG00000146648.1percent.500x.properpairs.2.bam"
    output_bam="/tmp/out.bam"
    main(input_bam, output_bam)
