import itertools
from pathlib import Path
import subprocess
from Bio.Data.IUPACData import ambiguous_dna_values, ambiguous_dna_complement
import seaborn as sns
import pandas as pd

def expand_ambiguous(seq):
    """
    >>> list(expand_ambiguous("NAW"))
    ['GAA', 'GAT', 'AAA', 'AAT', 'TAA', 'TAT', 'CAA', 'CAT']
    """
    for i in map("".join, itertools.product(*map(ambiguous_dna_values.get, seq))):
        yield i


def find_all_substrings(s, substr):
    """
    Find all (possibly overlapping) substrings in s and yield positions

    >>> list(find_all_substrings('banana', 'ana'))
    [1, 3]
    """
    i = s.find(substr)
    while i != -1:
        yield i
        i = s.find(substr, i + 1)


def bowtie_index(fasta, index_prefix):
    p = Path(index_prefix)
    p.mkdir(exist_ok=True, parents=True)
    cmds = ['bowtie-build', '-r', '-q', '-f', fasta, index_prefix]
    subprocess.run(cmds, check=True)


def bowtie(fasta, index_prefix, output, mismatches=3):
    cmds = [
        'bowtie',

        # report end-to-end hits with <= mismatches. Ignores qualities.
        '-v', str(mismatches),

        # report ALL hits
        '-a',

        # index
        '-x', index_prefix,

        # SAM output
        '-S',

        '--seed', '10', 

        # fasta, not fastq, input
        '-f', fasta
    ]
    with open(output, 'w') as fout:
        subprocess.run(cmds, stdout=fout)
    return output

def plot_coverage(annotation, grnas, outfn):
    """
    Plots the coverage by gRNAs for each feature in a single plot

    annotation : pybedtools object
        features
    grnas : pybedtools object
        candidates gRNAs
    outfn : str
        output name
    """
    # bedtools coverage per position
    cov = annotation.coverage(grnas, d=True).to_dataframe()
    cov.columns = ['contig', 'start', 'end', 'name', 'score', 'strand',
                    'feature_position', 'coverage']
    # for plotting the feature region
    cov['feature'] = 'gRNA'
    cov2 = cov.copy()
    cov2['feature'] = 'feature'
    cov2['coverage'] = -0.2
    cov = pd.concat([cov, cov2], ignore_index=True)

    # define figure parameters
    max_pos = max(cov['feature_position'])
    # plot each feature in annotation
    p = sns.FacetGrid(cov,
                      row="name",
                      hue='feature',
                      height=3,
                      aspect=max(max_pos/200, 2))
    p.map(sns.barplot, "feature_position", "coverage")
    # skip labels on x-axis for readability
    # choose a step depending on largest contig length
    stp = 10 if max_pos <= 100 else 50 if 100 < max_pos <= 500 else 100
    x_coord_range = range(stp-1, max_pos-1, stp)
    p.set(xticks=x_coord_range)
    p.set_xticklabels([x+1 for x in x_coord_range], rotation=90)
    p.add_legend()
    p.savefig(outfn + "_coverage.pdf")
