# 🧬 UVA SNV Pipeline

A SLURM-based whole-genome single-nucleotide variant (SNV) analysis pipeline for bacterial isolates. **Tested on UVA Rivanna HPC.**

## Overview

The pipeline performs:

- Read quality control and trimming
- Mash-based species identification and reference comparison for each isolate
- Two rounds of reference-guided SNV calling with Snippy
- SNP-based clustering with species annotations
- Reporting of SNP distance matrices, callable core-alignment metrics, and final group assignments
- Per-analysis software version and clustering-threshold provenance

### Why two rounds of variant calling?

Using one shared reference for a diverse collection of isolates can reduce the portion of the genome available for core-genome comparison and complicate interpretation of SNP distances. The pipeline therefore performs:

1. **Round 1:** SNV calling across all isolates using a selected shared reference.
2. **Initial clustering (default: 200 SNPs):** Identify broadly related groups.
3. **Round 2:** Select a reference for each group and repeat SNV calling within that group.
4. **Final clustering (default: 50 SNPs):** Resolve relatedness within each group.

Singleton/unrelated isolates are included in the final grouping table but are not reprocessed in Round 2. Comparing the callable alignment percentages between rounds helps reveal how divergent isolates affect the shared core genome.

## Installation

**Requirements:** A Linux HPC environment with SLURM, Conda (or a compatible environment manager), internet access for database installation, and storage accessible to both login and compute nodes. Commands below are run from the repository root unless stated otherwise.

### 1. Clone the repository and create the Conda environment

```bash
git clone https://github.com/scastanedabarba/uva-snvpipe.git
cd uva-snvpipe

# Create the environment on a filesystem accessible to compute nodes.
conda env create --prefix /path/to/shared/snvenv -f environment.yml
```

The `environment.yml` file defines the pipeline dependencies, including Snippy, Mash, read QC tools, SNP distance tools, and Python packages for alignment metrics. A site-specific environment prefix is recommended for shared HPC installations.

### 2. Configure the pipeline

Edit `config/config.sh` before installing the database or running an analysis:

```bash
# Absolute path to the installed Conda environment
CONDA_ENV="/path/to/shared/snvenv"

# Set after database installation (Step 3)
DB_DIR="/path/to/database_parent/databases/uva_eskape_2026-01-23"

# SLURM settings; leave empty to use cluster defaults
ACCOUNT="your_slurm_account"
PARTITION="your_partition"
```

Review the `CPUS_*`, `MEM_*`, and `TIME_*` settings in the same file to match your cluster's resources and limits. `config/config.sh` contains site-specific settings; analysis output locations are supplied with `--outdir` when launching the pipeline.

### 3. Download the reference database and set its path

The pipeline uses the **UVA ESKAPE reference database** (release `2026-01-23`).

- **Zenodo DOI:** [10.5281/zenodo.18838576](https://doi.org/10.5281/zenodo.18838576)
- **Archive:** `uva_eskape_2026-01-23.tar.gz`

The installer submits a SLURM job by default, downloads and verifies the archive, and extracts it into the `databases/` subdirectory of the supplied parent directory:

```bash
bash scripts/download_db.sh --outdir /path/to/database_parent
```

The resulting database path is:

```text
/path/to/database_parent/databases/uva_eskape_2026-01-23/
```

Check the installer's submitted job and logs for successful completion. Then set `DB_DIR` in `config/config.sh` to this **extracted database directory** (not the archive or its parent). The database contains `mash/`, `chrom_fna/`, and supporting metadata. The installer also supports `--local` for execution without SLURM submission.

### 4. Prepare a samplesheet and run the pipeline

Create a CSV file with the following header and one row per isolate:

```csv
sample,fastq_1,fastq_2
ISO1,/path/to/ISO1_R1.fastq.gz,/path/to/ISO1_R2.fastq.gz
ISO2,/path/to/ISO2_R1.fastq.gz,/path/to/ISO2_R2.fastq.gz
```

Use absolute paths to paired-end FASTQ files that compute nodes can access. Launch the pipeline with an explicit samplesheet and analysis output directory:

```bash
bash bin/run_pipeline.sh \
  --samplesheet /path/to/samplesheet.csv \
  --outdir /path/to/analysis_results
```

The launcher submits eight dependent SLURM stages. Follow their progress with `squeue -u "$USER"` or `sacct`, and inspect `slurm_logs/` in the analysis output directory if a stage fails. Successful job submission does not by itself confirm successful completion.

## Usage options

The database path defaults to `DB_DIR` in `config/config.sh`; override it for a particular run with `--db-path`:

```bash
bash bin/run_pipeline.sh \
  --samplesheet /path/to/samplesheet.csv \
  --outdir /path/to/analysis_results \
  --db-path /path/to/alternate_database
```

For a custom reference collection, you can specify the Mash sketch, mapping file, and reference FASTA directory explicitly:

```bash
bash bin/run_pipeline.sh \
  --samplesheet /path/to/samplesheet.csv \
  --outdir /path/to/analysis_results \
  --mash-sketch /path/to/database.msh \
  --mash-mapping /path/to/database.mapping.txt \
  --ref-dir /path/to/chrom_fna
```

The default clustering thresholds are 200 SNPs for Round 1 and 50 SNPs for Round 2. Override them if appropriate for the analysis:

```bash
bash bin/run_pipeline.sh \
  --samplesheet /path/to/samplesheet.csv \
  --outdir /path/to/analysis_results \
  --round1-threshold 200 \
  --round2-threshold 50
```

## Key outputs

All paths below are relative to the analysis directory passed to `--outdir`.

### Final group assignments

**`final_output/final_groups_all.tsv`** is the consolidated reporting table:

| Column | Description |
|---|---|
| `Isolate` | Sample identifier from the samplesheet |
| `Species` | Species annotation derived from Mash results |
| `Primary_Group` | Initial group assigned in Round 1 |
| `Secondary_Group` | Within-group assignment in Round 2 |

Unrelated/singleton isolates are retained and labeled `Unrelated` rather than reprocessed in Round 2. Group names are analysis-specific labels, not stable identifiers across different runs.

### SNP distance matrices

- **Round 1:** `final_output/round1/core.snp-dists.txt`
- **Round 2:** `final_output/round2/GroupX.core.snp-dists.txt` (one matrix per processed group)

These are pairwise SNP distance matrices from the respective core alignments. The **`Reference` row and column are intentionally retained** to show the selected reference's distance from each isolate. For relatedness interpretation, focus on isolate-to-isolate distances and consider the amount of shared callable genome.

The corresponding working outputs remain available under `variant_calling_round1/` and `variant_calling_round2/`.

### Core-alignment metrics

- **Round 1:** `final_output/round1/core_alignment_metrics.tsv`
- **Round 2:** `final_output/round2/GroupX.core_alignment_metrics.tsv`

Each TSV contains one summary row for its `core.full.aln`, with:

| Column | Meaning |
|---|---|
| `Alignment_Path` | Path to the source `core.full.aln` |
| `Total_Positions` | Total positions in the full alignment |
| `Callable_Positions` | Positions without `-`, `N`, `X`, or `n` in any non-reference isolate |
| `Excluded_Positions` | Positions excluded by that rule |
| `Percent_Callable` | Callable positions divided by total positions, as a percentage |

The percentage describes the shared callable alignment **across the isolates included in that round**, not an individual isolate's coverage. A divergent isolate can substantially reduce the Round 1 callable fraction; the Round 2 results provide a useful comparison after grouping more closely related isolates.

### Mash species identification and reference suitability

For each sample, review:

```text
mash/<sample>/<sample>.mash-screen.top3_hits.txt
```

This file reports the top Mash screen matches and should be checked for **species identification** and the strength of reference matches. As a practical QC guideline for this workflow, **more than 40,000 shared hashes** is considered supportive of an accurate reference mapping. This is a pipeline-specific heuristic rather than a universal Mash threshold: evaluate it alongside species concordance, competing matches, and the downstream alignment results. Low shared-hash counts or inconsistent species assignments warrant further review before interpreting SNP distances.

Additional Mash files are retained in each `mash/<sample>/` directory for troubleshooting and reference selection.

### Analysis provenance

**`analysis_provenance.tsv`** records the pipeline Git version/commit, installed software versions, database information, and the Round 1 and Round 2 clustering thresholds used for the analysis. It documents analysis configuration, not SLURM execution or completion status.

## Output directory structure

```text
analysis_results/
├── analysis_provenance.tsv
├── qc/
├── mash/
│   └── <sample>/
│       └── <sample>.mash-screen.top3_hits.txt
├── variant_calling_round1/
│   ├── findref/
│   ├── snippy/
│   └── summary/
├── variant_calling_round2/
│   ├── manifest/
│   ├── groups/
│   └── summary/
├── final_output/
│   ├── final_groups_all.tsv
│   ├── round1/
│   │   ├── core.snp-dists.txt
│   │   └── core_alignment_metrics.tsv
│   └── round2/
│       ├── Group1.core.snp-dists.txt
│       ├── Group1.core_alignment_metrics.tsv
│       └── ...
└── slurm_logs/
```

## Version

The pipeline's Git tag/commit and installed dependency versions are recorded per analysis in `analysis_provenance.tsv`. See repository releases for the current pipeline version.

