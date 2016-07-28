#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments,deprecated-method
from __future__ import print_function, absolute_import, division
import os
import pysam
from testfixtures.tempdirectory import TempDirectory

import test.utils_test as utils_test
from connor.samtools import BamFlag
import connor.samtools as samtools


def align(**kwargs):
    a = pysam.AlignedSegment()
    a.query_name = "align1"
    a.query_sequence="AGCTTAG"
    a.flag = 99
    a.reference_id = 0
    a.reference_start = 32
    a.mapping_quality = 20
    a.cigar = ((0,3), (2,1), (0,1))
    a.next_reference_id = 0
    a.next_reference_start=199
    a.template_length=167
    for (key, value) in kwargs.items():
        setattr(a, key, value)
    return a

def _create_file(path, filename, contents):
    filename = os.path.join(path, filename)
    with open(filename, 'wt') as new_file:
        new_file.write(contents)
        new_file.flush()
    return filename

def create_bam(path, filename, sam_contents, index=True):
    sam_filename = _create_file(path, filename, sam_contents)
    bam_filename = sam_filename.replace(".sam", ".bam")
    pysam_bam_from_sam(sam_filename, bam_filename, index)
    return bam_filename

def pysam_bam_from_sam(sam_filename, bam_filename, index=True):
    infile = samtools.alignment_file(sam_filename, "r")
    outfile = samtools.alignment_file(bam_filename, "wb", template=infile)
    for s in infile:
        outfile.write(s)
    infile.close()
    outfile.close()
    if index:
        samtools.index(bam_filename)

def pysam_alignments_from_bam(bam_filename):
    infile = samtools.alignment_file(bam_filename, "rb")
    aligned_segments = [s for s in infile]
    infile.close()
    return aligned_segments

class SamtoolsTest(utils_test.BaseConnorTestCase):
    def test_sort_and_index_bam(self):
        sam_contents = \
'''@HD|VN:1.4|GO:none|SO:coordinate
@SQ|SN:chr10|LN:135534747
readNameB1|147|chr10|400|0|5M|=|200|100|CCCCC|>>>>>
readNameA1|147|chr10|300|0|5M|=|100|100|AAAAA|>>>>>
readNameA1|99|chr10|100|0|5M|=|300|200|AAAAA|>>>>>
readNameB1|99|chr10|200|0|5M|=|400|200|CCCCC|>>>>>
readNameA2|147|chr10|300|0|5M|=|100|100|AAAAA|>>>>>
readNameA2|99|chr10|100|0|5M|=|300|200|AAAAA|>>>>>
'''.replace("|", "\t")

        with TempDirectory() as tmp_dir:
            bam = create_bam(tmp_dir.path,
                             "input.sam",
                             sam_contents,
                              index=False)
            samtools.sort_and_index_bam(bam)
            alignments = samtools.alignment_file(bam, "rb").fetch()
            aligns = [(a.query_name, a.reference_start + 1) for a in alignments]
            self.assertEquals(6, len(aligns))
            self.assertEquals([("readNameA1", 100),
                               ("readNameA2", 100),
                               ("readNameB1", 200),
                               ("readNameA1", 300),
                               ("readNameA2", 300),
                               ("readNameB1", 400)],
                              aligns)

            original_dir = os.getcwd()
            try:
                os.chdir(tmp_dir.path)
                os.mkdir("tmp")
                bam = create_bam(os.path.join(tmp_dir.path, "tmp"),
                                 "input.sam",
                                 sam_contents,
                                 index=False)
                bam_filename = os.path.basename(bam)

                samtools.sort_and_index_bam(os.path.join("tmp", bam_filename))

                alignments = samtools.alignment_file(bam, "rb").fetch()
                aligns = [(a.query_name,
                           a.reference_start + 1) for a in alignments]
                self.assertEquals(6, len(aligns))
                self.assertEquals([("readNameA1", 100),
                                   ("readNameA2", 100),
                                   ("readNameB1", 200),
                                   ("readNameA1", 300),
                                   ("readNameA2", 300),
                                   ("readNameB1", 400)],
                                  aligns)
            finally:
                os.chdir(original_dir)

    def test_filter_alignments_passthorughIncludedAligns(self):
        align1 = align(query_name="align1")
        base = [align1]

        filtered_aligns = [x for x in samtools.filter_alignments(base)]

        self.assertEqual(base, filtered_aligns)

    def test_filter_alignments_excludesUnpairedAligns(self):
        flag = 99
        align1 = align(query_name="align1", flag=flag)
        align2 = align(query_name="align2", flag=flag ^ BamFlag.PROPER_PAIR)
        align3 = align(query_name="align3", flag=flag)
        base = [align1, align2, align3]

        actual_names = [x.query_name for x in samtools.filter_alignments(base)]

        self.assertEqual(["align1", "align3"], actual_names)

    def test_filter_alignments_excludesSecondaryAligns(self):
        flag = 99
        align1 = align(query_name="align1", flag=flag)
        align2 = align(query_name="align2", flag=flag | BamFlag.SECONDARY)
        align3 = align(query_name="align3", flag=flag)
        base = [align1, align2, align3]

        actual_names = [x.query_name for x in samtools.filter_alignments(base)]

        self.assertEqual(["align1", "align3"], actual_names)

    def test_filter_alignments_excludesQCFails(self):
        flag = 99
        align1 = align(query_name="align1", flag=flag)
        align2 = align(query_name="align2", flag=flag | BamFlag.QCFAIL)
        align3 = align(query_name="align3", flag=flag)
        base = [align1, align2, align3]

        actual_names = [x.query_name for x in samtools.filter_alignments(base)]

        self.assertEqual(["align1", "align3"], actual_names)

    def test_filter_alignments_excludesMapq0(self):
        align1 = align(query_name="align1", mapping_quality=1)
        align2 = align(query_name="align2", mapping_quality=0)
        align3 = align(query_name="align3", mapping_quality=1)
        base = [align1, align2, align3]

        actual_names = [x.query_name for x in samtools.filter_alignments(base)]

        self.assertEqual(["align1", "align3"], actual_names)

    def test_filter_alignments_excludesCigarUnavailable(self):
        align1 = align(query_name="align1", cigarstring="6M")
        align2 = align(query_name="align2", cigarstring="*")
        align3 = align(query_name="align3", cigarstring="6M")
        base = [align1, align2, align3]

        actual_names = [x.query_name for x in samtools.filter_alignments(base)]

        self.assertEqual(["align1", "align3"], actual_names)

    def test_filter_alignments_logsFilterStats(self):
        align1 = align(query_name="align1", cigarstring="6M")
        align2 = align(query_name="align2", mapping_quality=0)
        align3 = align(query_name="align2", cigarstring="*", mapping_quality=0)
        align4 = align(query_name="align3", cigarstring="6M")
        base = [align1, align2, align3, align4]

        log = utils_test.MockLogger()
        for dummy in samtools.filter_alignments(base, log):
            pass

        self.assertEqual(1, len(log._log_calls['INFO']))
        self.assertEqual(log._log_calls['INFO'][0],
                        (r'2/4 (50.00%) alignments were excluded because they '
                         r'failed one or more filters (see log file for '
                         r'details)'))

        self.assertEqual(2, len(log._log_calls['DEBUG']))
        self.assertRegexpMatches(log._log_calls['DEBUG'][0],
                                 '1 alignment.*mapping quality')
        self.assertRegexpMatches(log._log_calls['DEBUG'][0],
                                 '1 alignment.*cigar unavail.*mapping quality')
