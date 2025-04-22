# Overview

This repo contains the source code corresponding to the full analysis in Luo et al., 2025.

The analyses are run through [Snakemake](https://snakemake.readthedocs.io/en/stable/)
pipelines. The corresponding python scripts, configuration files and downstream RMarkdown
files are included.

The pipelines are broken down into sections:

- design of DASH rRNA depletion oligos
- preparation of annotation files
- processing of RNAseq samples and Differential Analysis
- processing of RILseq samples and identification of RNA-RNA interactions
- processing of RIPseq samples





## Requirements

- [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) must
  be installed
- This code repository must be downloaded and unpacked. Below, we refer to the
  unpacked directory as `$WORKDIR`. Unless otherwise specified, paths provided
  are assumed to be relative to `$WORKDIR`.
- Either use a [snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
  appropriate for your cluster if running on a cluster, or run on a single machine (specifying the number
  of cores; here we use the placeholder $CORES).


## Create conda environment

Create the conda environments in the top level of `$WORKDIR`:

```bash
conda create -y -p env --file env-rilseq.yaml
```


## Run DASH rRNA depletion pipeline

This pipeline identifies depletion guide candidates (sgRNAs) targeted to specific genes, here
E. coli rRNA genes.

CITE THE OTHER - DASH - REPO. See OTHER REPO for usage and parameters.

Required inputs are a FASTA file, a GTF file of region(s) to deplete, a PAM sequence
(e.g. "NGG" for wildtype Cas9) and a guide length. The `workflows/ecoli-dash/config/config.yaml`
file configures this. `workflows/ecoli-dash/config/rrna_names.tsv` contains the rRNA gene names
patterns. Run the Snakefile pipeline to download the references, format them for DASH and run
`CandidateFinder`.

In the `workflows/ecoli-dash` directory:

```bash
# activate the environment
conda activate ../../env

# run the DASH pipeline
snakemake --cores=$CORES --use-conda
```

This generates the output file `workflows/ecoli-dash/data/ecoli.GCF_000005845.2_ASM584v2_20/candidates_oligos.csv`,
a comma-separated file of unique full oligo sequences. They are built from the T7 promoter
sequence, the 20-bp candidate guides, the Cas9 sgRNA scaffold and the reverse complement of
the primer.


## Run annotation pipeline

The annotation Snakemake pipeline will download the references, generate the aligners indexes
and reformat the annotation to contain the CDS features, 5UTR, 3UTR, AS (antisense) and IGR
(intergenic regions). The `workflows/annotation/config/config.yaml` file configures this.

To run the workflow, navigate to the `workflows/annotation` directory, then activate the
environment and run the Snakefile with the following commands:

```bash
# activate the environment
conda activate ../../env

# run the annotation pipeline
snakemake --cores $CORES
```

The main outputs are the RILseq-formatted annotation file, `workflows/annotation/data/ecoli_default.patch3.gff`,
where each feature (CDS, UTRs, AS, IGR) are assigned their corresponding feature, and the
RNAseq-formatted annotation file, `workflows/annotation/data/ecoli_default.patch4.gff` where
every feature is assigned "exon" so they can be quantified by featureCounts. These are already configured
to be used in other workflows, in their respective `workflows/{workflow}/config/config.yaml` files.

???? and the bwa index file, `workflows/annotation/data/ecoli_default.fasta.bwt`. 


## Run RNAseq pipeline

All RNAseq samples are processed through a single RNAseq pipeline. Downstream analyses are
split in several downstream RMarkdown files. The `workflows/rnaseq-FMD-AX/config/config.yaml`
file has already been configured. It specifies parameters such as the genome reference to download,
the previously generated annotation file, which aligner to use. See
[lcdb-wf configuration](https://lcdb.github.io/lcdb-wf/config.html) for details.

To run the workflow, navigate to the `workflows/rnaseq-FMD-AX` directory, then activate the
environment and run the Snakefile with the following commands:

```bash
# activate the environment
conda activate ../../env

# run the annotation pipeline
snakemake --cores $CORES --use-conda
```






