#pylint: disable=invalid-name, too-few-public-methods, too-many-public-methods
#pylint: disable=protected-access, missing-docstring, too-many-locals
#pylint: disable=too-many-arguments,deprecated-method
from __future__ import print_function, absolute_import, division
import os
from testfixtures.tempdirectory import TempDirectory

from test.utils_test import BaseConnorTestCase
import connor.samtools as samtools

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

class SamtoolsTest(BaseConnorTestCase):
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
