'''
Created on Jun 3, 2016

@author: pulintz, cgates
'''
from __future__ import print_function, absolute_import, division
from _collections import defaultdict
from Bio.AlignIO.Interfaces import AlignmentWriter

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
    start_alignment = alignments[0]
    for alignment in alignments[1:]:
        if alignment.query_name == start_alignment.query_name:
            return (start_alignment, alignment)


if __name__ == '__main__':
    pass
