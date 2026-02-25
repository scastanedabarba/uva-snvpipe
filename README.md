# UVA SNV Pipeline

A SLURM-based whole genome SNV analysis pipeline for bacterial isolates (tested on UVA Rivanna).

This pipeline performs:

- Read QC and trimming
- Mash-based speciation (per isolate)
- Two rounds of reference-guided SNV calling (Snippy)
- SNP-based clustering with species annotation

## Why two rounds of variant calling?

Using a **single, shared reference** for a mixed set of isolates can inflate SNP distances when one or more isolates are very distant (mapping/variant calling becomes less comparable across the cohort).

To reduce this effect, the pipeline runs:

1. **Variant calling Round 1**: call SNPs across all isolates and compute a distance matrix.
2. **Initial clustering (threshold = 200 SNP)**: define *broad* relatedness groups.
3. **Variant calling Round 2 (within each broad group)**: re-select a reference *per group*, re-run Snippy and snippy-core, then compute **final** SNP distances and clusters within each group. The default threshold for final SNP clustering can be specified by the user, 50 is set as the default.

Singletons / unrelated isolates (no links under the Round 1 threshold) are carried through as **Unrelated** and are not re-run in Round 2.

---

# Installation

## 1. Clone the repository

```bash
git clone https://github.com/scastanedabarba/uva-snvpipe.git
cd uva-snvpipe
```

## 2. Configure environment

This pipeline runs using HPC modules and a conda environment.

Required tools:

- Mash 2.3
- TrimGalore
- Trimmomatic 0.39
- Java
- Snippy
- snp-sites
- snp-dists
- Python 3

Edit:

```bash
config/config.sh
```

Important variables:

```bash
ACCOUNT=
PARTITION=
OUTDIR=
MASH_DB=/path/to/pathogen_refseq.chrom.k21s50000.msh
MODE=modules
```

> Note: This version is **modules-only** (no containers).

---

# Samplesheet Format

CSV header:

```
sample,fastq_1,fastq_2
```

Example:

```
PGCoE1,/scratch/reads/PGCoE1_R1.fastq.gz,/scratch/reads/PGCoE1_R2.fastq.gz
PGCoE2,/scratch/reads/PGCoE2_R1.fastq.gz,/scratch/reads/PGCoE2_R2.fastq.gz
```

Absolute paths to FASTQ files are required.

---

# Usage

## Run with defaults from config.sh

```bash
bash bin/run_pipeline.sh
```

## Specify a samplesheet

```bash
bash bin/run_pipeline.sh --samplesheet config/samplesheet.csv
```

## Specify an output directory

```bash
bash bin/run_pipeline.sh --samplesheet config/samplesheet.csv --outdir /scratch/sgj4qr/snv_pipeline/test_run
```

If `--outdir` is provided, it overrides `OUTDIR` in `config.sh`.


## Specify thresholds for clustering in each round of varint calling (step 5 and 8)

```bash
bash bin/run_pipeline.sh --samplesheet config/samplesheet.csv --outdir /scratch/sgj4qr/snv_pipeline/test_run 
  --round1-threshold 200 
  --round2-threshold 50
```

---

# Pipeline Workflow (high level)

## Step 1 — QC

- TrimGalore
- Trimmomatic
- Produces trimmed paired FASTQs

Output:

```
OUTDIR/qc/<sample>/
```

## Step 2 — Mash speciation (per isolate)

- Concatenates trimmed reads
- `mash sketch`, `mash dist`, `mash screen`
- Writes top hits

Output:

```
OUTDIR/mash/<sample>/
```

### Where speciation results are stored

Primary speciation source (top hit):

```
OUTDIR/mash/<sample>/<sample>.mash-screen.top3_hits.txt
```

Fallback:

```
OUTDIR/mash/<sample>/<sample>.mash-ref.txt
```

The clustering script reduces the Mash description to **Genus species** (first two words).

---

## Variant calling Round 1 (all isolates)

### Step 3 — Reference selection (Round 1)

Selects a cohort reference that minimizes total Mash triangle distance across isolates.

Output:

```
OUTDIR/variant_calling_round1/findref/ref.fa
```

### Step 4 — SNV calling (Round 1)

Runs Snippy per isolate using the Round 1 reference.

Output:

```
OUTDIR/variant_calling_round1/snippy/<sample>/
```

### Step 5 — Core genome + initial clustering (Round 1)

Runs:

- `snippy-core`
- `snp-sites`
- `snp-dists`
- **Initial clustering using threshold = 200 SNP**

Core artifacts (Round 1):

```
OUTDIR/variant_calling_round1/snippy/core.aln
OUTDIR/variant_calling_round1/snippy/core.full.aln
OUTDIR/variant_calling_round1/snippy/phylo.aln
OUTDIR/variant_calling_round1/snippy/snippy_core.log
```

Round 1 SNP matrix + initial groups:

```
OUTDIR/variant_calling_round1/summary/core.snp-dists.txt
OUTDIR/variant_calling_round1/summary/initial_groups.tsv
```

---

## Variant calling Round 2 (within initial groups)

Round 2 only runs for initial groups with **≥ 2 isolates**. Singletons/unrelated isolates are not reprocessed.

### Step 6 — Reference selection (Round 2, per group)

Creates a Round 2 manifest and chooses a reference **per initial group**.

Manifest:

```
OUTDIR/variant_calling_round2/manifest/groups.tsv
OUTDIR/variant_calling_round2/manifest/group_members/GroupX.samples.txt
OUTDIR/variant_calling_round2/manifest/singletons.tsv
```

Per-group reference:

```
OUTDIR/variant_calling_round2/groups/GroupX/findref/ref.fa
```

### Step 7 — SNV calling (Round 2, per group)

Runs Snippy for group members using that group’s reference.

Per-group outputs:

```
OUTDIR/variant_calling_round2/groups/GroupX/snippy/<sample>/
```

### Step 8 — Core genome + final clustering (Round 2, per group)

Runs snippy-core per group and clusters **within-group** at **threshold = 50 SNP**.

Core artifacts (Round 2 are stored *within the group snippy directory* to match Round 1 layout):

```
OUTDIR/variant_calling_round2/groups/GroupX/snippy/core.aln
OUTDIR/variant_calling_round2/groups/GroupX/snippy/core.full.aln
OUTDIR/variant_calling_round2/groups/GroupX/snippy/phylo.aln
OUTDIR/variant_calling_round2/groups/GroupX/snippy/snippy_core.log
```

Round 2 SNP matrices + final group tables are written centrally:

```
OUTDIR/variant_calling_round2/summary/GroupX.core.snp-dists.txt
OUTDIR/variant_calling_round2/summary/GroupX.final_groups.tsv
```

---

# Key Outputs to use

## 1) Final table (recommended starting point)

**This is the main deliverable**: isolates with species + final grouping.

```
OUTDIR/final_output/final_groups_all.tsv
```

It includes:

- `Isolate`
- `Species` (Genus species)
- `Primary_Group` (Round 1 broad group; threshold 200)
- `Secondary_Group` (Round 2 within-group cluster; threshold 50; `Unrelated` if not reprocessed)

## 2) SNP distance matrices (variant matrices)

- Round 1 matrix (all isolates):

```
OUTDIR/variant_calling_round1/summary/core.snp-dists.txt
```

- Round 2 matrices (per initial group):

```
OUTDIR/variant_calling_round2/summary/GroupX.core.snp-dists.txt
```

The pipeline also copies Round 2 matrices into:

```
OUTDIR/final_output/GroupX.core.snp-dists.txt
```

---

# Output Directory Structure

```
OUTDIR/
├── qc/
├── mash/
├── variant_calling_round1/
│   ├── findref/
│   ├── snippy/
│   │   ├── core.aln
│   │   ├── core.full.aln
│   │   ├── phylo.aln
│   │   ├── snippy_core.log
│   │   └── <sample>/
│   ├── summary/
│   │   ├── core.snp-dists.txt
│   │   └── initial_groups.tsv
│   └── status/
│
├── variant_calling_round2/
│   ├── manifest/
│   │   ├── groups.tsv
│   │   ├── group_members/
│   │   └── singletons.tsv
│   ├── groups/
│   │   └── GroupX/
│   │       ├── samplesheet.csv
│   │       ├── findref/
│   │       └── snippy/
│   │           ├── core.aln
│   │           ├── core.full.aln
│   │           ├── phylo.aln
│   │           ├── snippy_core.log
│   │           └── <sample>/
│   └── summary/
│       ├── Group1.core.snp-dists.txt
│       ├── Group1.final_groups.tsv
│       ├── Group2.core.snp-dists.txt
│       ├── Group2.final_groups.tsv
│       └── ...
│
├── final_output/
│   ├── final_groups_all.tsv
│   ├── Group1.core.snp-dists.txt
│   ├── Group2.core.snp-dists.txt
│   └── ...
│
└── slurm_logs/
```

---

# Version

v0.2.0  
Two-round SNV calling + grouping (200 SNP initial, 50 SNP final within-group), with user-friendly `final_output/`.

