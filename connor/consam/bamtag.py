'''Defines tags to be added to output BAM file(s)'''
from __future__ import print_function, absolute_import, division

class BamTag(object):
    HEADER_FORMAT = 'connor|BAM tag|{}: {}'.replace('|', '\t')

    class _NullObject(object):
        '''Returns None for all method calls'''
        #pylint: disable=no-self-use, unused-argument
        def __init__(self):
            self.included_pair_count = None
            self.filter_value = None
            self.umi_sequence = None
            self.umt=lambda *args: None
            self.is_consensus_template = lambda *args: None
            self.positions = lambda *args: None
            self.cigars = lambda *args: None

    _NULL_OBJECT = _NullObject()

    def __init__(self, tag_name, tag_type, description, get_value):
        self._tag_name = tag_name
        self._tag_type = tag_type
        self._get_value = get_value
        self._description = description
        self.header_comment = BamTag.HEADER_FORMAT.format(tag_name,
                                                          description)

    def __lt__(self, other):
        return (self._tag_name,
                self._description) < (other._tag_name, other._description)

    def set_tag(self, family, paired_align, connor_align):
        family = family if family else BamTag._NULL_OBJECT
        paired_align = paired_align if paired_align else BamTag._NULL_OBJECT
        value = self._get_value(family, paired_align, connor_align)
        connor_align.set_tag(self._tag_name, value, self._tag_type)


def build_bam_tags():
    #pylint: disable=unused-argument
    def combine_filters(fam, paired_align, align):
        filters = [x.filter_value for x in [fam, align] if x and x.filter_value]
        if filters:
            return ";".join(filters).replace('; ', ';')
        else:
            return None
    boolean_tag_value = {True:1}
    tags = [
        BamTag("X0", "Z",
               ("filter (why the alignment was excluded)"),
               combine_filters),
        BamTag("X1", "Z",
               ("leftmost~rightmost matched pair positions"),
               lambda fam, pair, align: pair.positions('{left}~{right}')),
        BamTag("X2", "Z",
               ("L~R CIGARs"),
               lambda fam, pair, align: pair.cigars('{left}~{right}')),
        BamTag("X3", "i",
               "unique identifier for this alignment family",
               lambda fam, pair, align: fam.umi_sequence),
        BamTag("X4", "Z",
               ("L~R UMT barcodes for this alignment family; because "
                "of fuzzy matching the family UMT may be distinct "
                "from the UMT of the original alignment"),
               lambda fam, pair, align: fam.umt('{left}~{right}')),
        BamTag("X5", "i",
               "family size (number of align pairs in this family)",
               lambda fam, pair, align: fam.included_pair_count),
        BamTag("X6", "i",
               ("presence of this tag signals that this alignment "
                "would be the template for the consensus alignment"),
               lambda fam, pair, align: boolean_tag_value.get(fam.is_consensus_template(align), None))]
    return tags
