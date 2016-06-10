'''
Created on Jun 3, 2016

@author: pulintz, cgates
'''
from __future__ import print_function, absolute_import, division
from collections import defaultdict
import sys
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

class PairedAlignment(object):
    def __init__(self, left_alignment, right_alignment):
        self.left_alignment = left_alignment
        self.right_alignment = right_alignment

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):
        return hash(self.left_alignment) * hash(self.right_alignment)

    def __repr__(self):
        return ("Pair({}|{}|{}, "
                "{}|{}|{})").format(self.left_alignment.query_name, 
                                    self.left_alignment.reference_start,
                                    self.left_alignment.sequence,
                                    self.right_alignment.query_name, 
                                    self.right_alignment.reference_start,
                                    self.right_alignment.sequence)

def _build_coordinate_read_name_manifest(lw_aligns):
    af_dict = defaultdict(set)
    for lwa in lw_aligns:
        af_dict[lwa.key].add(lwa.name)
    return af_dict

def _build_coordinate_families(aligned_segments,coord_read_name_manifest):
    family_dict = defaultdict(set)
    pairing_dict = {}
    for aseg in aligned_segments:
        if not aseg.query_name in pairing_dict:
            pairing_dict[aseg.query_name]= aseg
        else:
            paired_align = PairedAlignment(pairing_dict.pop(aseg.query_name),
                                           aseg)
            key = LightweightAlignment(aseg).key
            family_dict[key].add(paired_align)
            coord_read_name_manifest[key].remove(aseg.query_name)
            if not coord_read_name_manifest[key]:
                yield family_dict.pop(key)

def _build_consensus_pair(alignment_family):
    return alignment_family.pop()


def _get_ranked_tags(coordinate_family):
    tag_count = defaultdict(int)
    for paired_align in coordinate_family:
        left_tag_id = paired_align.left_alignment.sequence[0:3]
        right_tag_id = paired_align.right_alignment.sequence[0:3]
        tag_count[(left_tag_id, right_tag_id)] += 1
    tags_by_count = sorted(tag_count.items(),
                           key=lambda x: (x[1], x[0]),
                           reverse=True)
    ranked_tags = [tag_count[0] for (tag_count) in tags_by_count]
    return ranked_tags


def _build_tag_families(coordinate_family):
    ranked_tags = _get_ranked_tags(coordinate_family)
    tag_aligns = defaultdict(set)
    for paired_align in coordinate_family:
        left_tag_id =  paired_align.left_alignment.sequence[0:3]
        right_tag_id =  paired_align.right_alignment.sequence[0:3]
        for best_tag in ranked_tags:
            if left_tag_id == best_tag[0] or right_tag_id == best_tag[1]:
                tag_aligns[best_tag].add(paired_align)
                break
    return tag_aligns.values()


def main(input_bam, output_bam):
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    lw_aligns = [LightweightAlignment(align) for align in bamfile.fetch()]
    coord_manifest = _build_coordinate_read_name_manifest(lw_aligns)
    bamfile.close()
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    outfile = pysam.AlignmentFile(output_bam, "wb", template=bamfile)
    for family in _build_coordinate_families(bamfile.fetch(), coord_manifest):
        read_pair = _build_consensus_pair(family)
        outfile.write(read_pair.left_alignment)
        outfile.write(read_pair.right_alignment)
    outfile.close()
    bamfile.close()

if __name__ == '__main__':
    in_bam = sys.argv[1]
    out_bam = sys.argv[2]
    main(in_bam, out_bam)
