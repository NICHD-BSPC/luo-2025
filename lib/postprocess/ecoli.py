from Bio import SeqIO
from snakemake import shell
import gzip
import os

def to_chr(origfn, newfn):
    """
    Asserts that there is only one sequence in the FASTA and renames it to
    "chr".
    """
    assert len(origfn) == 1

    tmp = origfn[0]+'.tmp.gz'
    shell('gunzip -c {origfn} > {tmp}')
    recs = list(SeqIO.parse(tmp, 'fasta'))
    assert len(recs) == 1
    recs[0].name = 'chr'
    recs[0].id = 'chr'
    recs[0].description = 'chr'
    with gzip.open(newfn, 'wt') as fout:
        SeqIO.write(recs, fout, 'fasta')
    os.unlink(tmp)
