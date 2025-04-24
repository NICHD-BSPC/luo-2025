from Bio import Seq
#
# T7 promoter (minus the first G)
T7 = "TTCTAATACGACTCACTATA"

# Cas9 sgRNA scaffold
SCAFFOLD = (
    "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"
)

# Reverse complement of this gets added to template sgRNA oligos
PRIMER = "GGCATACTCTGCGACATCGT"
PRIMER_REVCOMP = str(Seq.Seq(PRIMER).reverse_complement())

# Match output to original
LEGACY = False
