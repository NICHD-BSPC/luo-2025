# Overview

This repo contains the source code corresponding to the full analysis in Luo et al., 2025.

The analyses are run through [Snakemake](https://snakemake.readthedocs.io/en/stable/)
pipelines. The corresponding python scripts, configuration files and downstream RMarkdown
files are included.

The pipelines are broken down into sections:

- design of DASH rRNA depletion oligos
- preparation of annotation files
- processing of RNAseq samples, Differential Expression Analysis of Total RNA samples and
  RIP analysis of IP and corresponding input samples
- processing of RILseq samples and identification of RNA-RNA interactions


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
conda env create -y -p env --file env.yaml
```


## Run DASH rRNA depletion workflow

This workflow identifies depletion guide candidates (sgRNAs) targeted to specific genes, here
E. coli rRNA genes. It is based on the Github repository
[bacteria-dash](https://github.com/NICHD-BSPC/bacteria-dash).
A Snakefile downloads the references, formats them for DASH processing and runs `CandidateFinder`.
This will identify all possible PAM sequences on either strands of the regions to deplete - here rRNAs -
and returns candidates sgRNA of the specified length that have no off-target match.

Required inputs are the genome FASTA file, the annotation GTF file, a tsv of gene name pattern(s) to deplete,
the PAM sequence (e.g. "NGG" for wildtype Cas9) and the guide length to generate. The
`workflows/ecoli-dash/config/config.yaml` file configures this. `workflows/ecoli-dash/config/rrna_names.tsv`
contains the rRNA gene names patterns.

To run the Snakefile pipeline, navigate to the `workflows/ecoli-dash` directory, then activate the
environment and run the Snakefile with the following commands:

```bash
# activate the environment
conda activate ../../env

# run the DASH workflow
snakemake --cores=$CORES --use-conda
```

This generates the output file `workflows/ecoli-dash/data/ecoli.GCF_000005845.2_ASM584v2_20/candidates_oligos_annotated.tsv`,
a comma-separated file of unique full oligo sequences. They are built from the T7 promoter
sequence, the 20-bp candidate guides, the Cas9 sgRNA scaffold and the reverse complement of
the primer. Targeted rRNAs are indicated as well.


## Run annotation workflow

This workflow will download the references, generate the BWA indexes and reformat the annotation
to contain the CDS features, 5UTR, 3UTR, AS (antisense) and IGR (intergenic regions) required
for the RILseq and RNAseq workflows. The `workflows/annotation/config/config.yaml`
file configures this.

To run the workflow, navigate to the `workflows/annotation` directory, then activate the
environment and run the Snakefile with the following commands:

```bash
# activate the environment
conda activate ../../env

# run the annotation workflow
snakemake --cores $CORES
```

The main outputs are the RILseq-formatted annotation file, `workflows/annotation/data/ecoli_default.patch3.gff`,
where all individual features (CDS, UTRs, AS, IGR) are retained as their corresponding feature, the
RNAseq-formatted annotation file, `workflows/annotation/data/ecoli_default.patch4.gff` where
every feature is assigned "exon" so they can be quantified by featureCounts. This also creates
the bwa index file, `workflows/annotation/data/ecoli_default.fasta.bwt`. Paths to those files are
already configured to be used in other workflows, in their respective `workflows/{workflow}/config/config.yaml`
files.


## Run RNAseq workflow

Prerequisite: annotation workflow should be run first to produce the annotation file.

All RNAseq samples are processed through a single snakemake pipeline, based on the
[lcdb-wf](https://github.com/lcdb/lcdb-wf) github repository. The initial steps of QC filtering
and alignment to the genome are common to Total RNA (RNAseq) and IP or input (RIPseq) samples.

Downstream analyses are split by RNAseq Differential Expression and RIP differential analysis, and
are run through the Snakefile. RNAseq Differential Expression contrasts are performed between Total
RNA samples of different conditions; RIP analysis contrasts are performed between IP samples vs.
input samples.

The `workflows/rnaseq/config/config.yaml` and `workflows/rnaseq/config/sampletable.tsv`
files have already been configured, but can be edited to match other samples and experimental design.
The config.yaml file specifies parameters such as the genome reference to download, the annotation file,
which aligner to use. See [lcdb-wf configuration](https://lcdb.github.io/lcdb-wf/config.html)
for details.

Here we use the annotation generated above by the annotation workflow, which contains in addition
to exons the UTR, AS, and IGR features, so every feature type can be quantified in `featureCounts`.
Note that the `featureCounts` parameter `-O` to allow counting of reads overlapping multiple features
was required to include reads spanning UTRs and exons.

To run the workflow, navigate to the `workflows/rnaseq` directory, then activate the
environment and run the Snakefile with the following commands:

```bash
# activate the environment
conda activate ../../env

# run the rnaseq workflow
snakemake --cores $CORES --use-conda
```

Workflow outputs include the raw read counts, which can be found in
`workflows/rnaseq/data/rnaseq_aggregation/featurecounts_all_features_allow_overlap.txt`.

The RNAseq differential expression and RIP differential analyses main outputs are the individual  html
reports, for each of the analysis subsets, `workflows/rnaseq/rnaseq-downstream-{subset}/rnaseq.html`
and `workflows/rnaseq/ripseq-downstream-{subset}/rnaseq.html`. In those reports are PCA plots,
volcano plots, and tables of differentially expressed features. Links to high resolution PDFs and
TSV tables are provided.


## Run RILseq workflow

Prerequisite: the annotation workflow and RNAseq workflows should be run first to produce the annotation
and counts files respectively.

This workflow identifies RNA-RNA interactions by finding chimera transcripts in cross-linked samples.
The snakemake pipeline starts by generating the cutadapt-trimmed FASTQ files using `cutadapt`,
then identify single and chimeric fragments among the reads and filter them for significant chimeric fragments.
It uses functions from the [RIL-seq package](https://github.com/asafpr/RILseq). The Snakefile also
runs downstream Rmarkdown scripts that generate summary figures and tables of the mapped fragments. Note that
this step frequently takes several hours.

The sample table (`workflows/rilseq/config/sampletable.tsv`) and configuration file
(`workflows/rilseq/config/config.yaml`) are already configured, but can be edited to match other samples
and experimental design.

Circos plots, which show the significant chimera locations, need a configuration file that matches
the reference genome and annotation. Here we used the default E. coli file provided in the
[RIL-seq](https://github.com/asafpr/RILseq) Github repo, which can be found in
`workflows/rilseq/config/RILseq-master/data/E_coli_K12/ver2/draw_links.conf`.

To run the workflow, navigate to the `workflows/rilseq` directory, , then activate the
environment and run the Snakefile with the following commands:

```bash
# activate the environment
conda activate ../../env

# run the rilseq workflow
snakemake --cores $CORES --use-conda
```

The summary reports and xlsx tables aggregating the results are split by groups of samples, for
convenience and to take advantage of the snakemake parallelisation. The main outputs are
the summary reports in `workflows/rilseq/downstream-{subset}/summary.html`. They contain various
QC and interaction quantification plots, and links to the xlsx aggregation tables
`workflows/rilseq/downstream-{subset}/Summary_RILSeq_S-chimeras_{subset}.xlsx` and
`workflows/rilseq/downstream-{subset}/Summary_interactions_with_IP_{subset}.xlsx`.

`Summary_RILSeq_S-chimeras_{subset}.xlsx` contains for each individual sample the positions and
counts of chimeras supporting significant interactions, and their corresponding statistical test results.
Additional spreadsheets contain the interactins found in common between replicates.

`Summary_interactions_with_IP_{subset}.xlsx` contains the chimera counts per features for the chimeric
reads, single-transcript reads, and when applicable the reads in their corresponding Total RNA samples.

Replicates can be aggregated by merging the reads and processing the resulting BAM files as their
individual replicate counterparts. They are designated as `unified` samples. The results can be
found in `workflows/rilseq/downstream-{subset}/summary-unified.html` and 
`workflows/rilseq/downstream-{subset}/Summary_RILSeq_S-chimeras_unified_{subset}.xlsx`. 


