# UVA SNV Pipeline

A SLURM-based whole genome SNV analysis pipeline for bacterial isolates.

This pipeline performs:

1.  Read QC and trimming\
2.  Mash-based speciation and reference selection\
3.  Reference-guided SNV calling (Snippy)\
4.  Core genome alignment\
5.  SNP distance matrix generation\
6.  SNP-based clustering (≤ 50 SNP threshold) with species annotation

Designed for HPC environments (tested on UVA Rivanna).

------------------------------------------------------------------------

# Installation

## 1. Clone the repository

``` bash
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>
```

## 2. Configure environment

This pipeline currently runs using **HPC modules and conda
environments**.

Ensure the following tools are available:

-   Mash 2.3
-   TrimGalore
-   Trimmomatic 0.39
-   Java
-   Snippy
-   snp-sites
-   snp-dists
-   Python 3

Example (Rivanna-style):

``` bash
module use /project/amr_services/modulefiles/
module load mash/2.3
module load trimgalore
module load trimmomatic/0.39
module load snp-dists
module load miniforge
source activate /project/amr_services/.conda/amrpipe-snippy-env
```

Modify paths in:

``` bash
config/config.sh
```

Important variables:

``` bash
ACCOUNT=
PARTITION=
OUTDIR=
MASH_DB=/path/to/pathogen_refseq.chrom.k21s50000.msh
```

------------------------------------------------------------------------

# Samplesheet Format

CSV file with header:

    sample,fastq_1,fastq_2

Example:

    PGCoE1,/scratch/reads/PGCoE1_R1.fastq.gz,/scratch/reads/PGCoE1_R2.fastq.gz
    PGCoE2,/scratch/reads/PGCoE2_R1.fastq.gz,/scratch/reads/PGCoE2_R2.fastq.gz

Absolute paths to FASTQ files are recommended.

------------------------------------------------------------------------

# Usage

## Basic Run (using defaults in config.sh)

``` bash
bash bin/run_pipeline.sh
```

## Specify Samplesheet

``` bash
bash bin/run_pipeline.sh   --samplesheet config/samplesheet.csv
```

## Specify Output Directory

``` bash
bash bin/run_pipeline.sh   --samplesheet config/samplesheet.csv   --outdir /scratch/sgj4qr/snv_pipeline/test_run
```

If `--outdir` is provided, it overrides `OUTDIR` in `config.sh`.

------------------------------------------------------------------------

# Pipeline Workflow

## Step 1 -- QC

-   TrimGalore
-   Trimmomatic
-   Generates trimmed paired FASTQs

Output:

    OUTDIR/qc/<sample>/

------------------------------------------------------------------------

## Step 2 -- Mash Speciation

-   Concatenates trimmed reads
-   mash sketch
-   mash dist
-   mash screen
-   Extracts top hit

Output:

    OUTDIR/mash/<sample>/

------------------------------------------------------------------------

## Step 3 -- Reference Selection

-   Copies isolate sketches
-   Sketches candidate reference genomes
-   Runs mash triangle
-   Selects reference minimizing total distance

Writes:

    OUTDIR/findref/ref.fa

------------------------------------------------------------------------

## Step 4 -- SNV Calling

-   Runs Snippy per isolate
-   Uses selected reference

Output:

    OUTDIR/snippy/<sample>/

------------------------------------------------------------------------

## Step 5 -- Core Genome + Clustering

Runs:

-   snippy-core
-   snp-sites
-   snp-dists
-   SNP-based clustering (≤ 50 SNP)

Outputs:

    OUTDIR/snippy/core.aln
    OUTDIR/snippy/core.full.aln
    OUTDIR/snippy/phylo.aln
    OUTDIR/snippy/snippy_core.log

    OUTDIR/summary/core.snp-dists.txt
    OUTDIR/summary/isolate_groups.tsv

------------------------------------------------------------------------

# Clustering Logic

-   Graph-based clustering
-   Isolates connected if SNP distance ≤ 50
-   Connected components become:

```{=html}
<!-- -->
```
    Group1
    Group2
    ...
    Unrelated

Species is extracted from the first Mash screen hit and reduced to:

    Genus species

Final output file:

    summary/isolate_groups.tsv

Columns:

    Isolate    Species    Group

------------------------------------------------------------------------

# Output Directory Structure

    OUTDIR/
    ├── qc/
    ├── mash/
    ├── findref/
    ├── snippy/
    │   ├── core.aln
    │   ├── core.full.aln
    │   ├── phylo.aln
    │   ├── snippy_core.log
    │   └── <per isolate directories>
    │
    ├── summary/
    │   ├── core.snp-dists.txt
    │   └── isolate_groups.tsv
    │
    ├── status/
    └── slurm_logs/

------------------------------------------------------------------------

# Current Version

v0.1.0\
First fully operational SLURM + modules-based release.

