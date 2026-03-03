# 🧬 UVA SNV Pipeline

A SLURM-based whole genome SNV analysis pipeline for bacterial isolates\
(Tested on UVA Rivanna HPC)

------------------------------------------------------------------------

## 🚀 Overview

This pipeline performs:

-   Read QC and trimming
-   Mash-based speciation
-   **Two rounds of reference-guided SNV calling (Snippy)**
-   SNP-based clustering with species annotation
-   Final consolidated grouping table for reporting

------------------------------------------------------------------------

```{=html}
<details>
```
```{=html}
<summary>
```
`<strong>`{=html}🧠 Why Two Rounds of Variant Calling?`</strong>`{=html}
```{=html}
</summary>
```
Using a single shared reference for a diverse set of isolates can
inflate SNP distances when distant isolates are included.

To improve clustering resolution:

1.  **Round 1** --- SNV calling across all isolates\
2.  **Initial clustering (default: 200 SNP)** --- define broad
    relatedness groups\
3.  **Round 2** --- Re-select a reference per group and re-run SNV
    calling\
4.  **Final clustering (default: 50 SNP)** --- within-group resolution

Singleton/unrelated isolates are carried forward but not reprocessed in
Round 2.

```{=html}
</details>
```

------------------------------------------------------------------------

# 📦 Installation

## 1️⃣ Clone the repository

``` bash
git clone https://github.com/scastanedabarba/uva-snvpipe.git
cd uva-snvpipe
```

------------------------------------------------------------------------

## 2️⃣ Download the Reference Database (Required)

The pipeline requires the **UVA ESKAPE reference database**.

**Zenodo DOI:**\
`10.5281/zenodo.18838576`

### Install to default location (recommended)

``` bash
bash scripts/download_db.sh
```

Installs to:

    databases/uva_eskape_2026-01-23/

### Install to a custom location

``` bash
bash scripts/download_db.sh --output /scratch/my_databases
```

Then run with:

``` bash
bash bin/run_pipeline.sh --db-path /scratch/my_databases/uva_eskape_2026-01-23
```

------------------------------------------------------------------------

## 3️⃣ Configure Environment

Edit:

``` bash
config/config.sh
```

Set:

``` bash
ACCOUNT=
PARTITION=
OUTDIR=
MODE=modules
```

### Required Tools

-   Mash 2.3
-   TrimGalore
-   Trimmomatic 0.39
-   Java
-   Snippy
-   snp-sites
-   snp-dists
-   Python 3

> ⚠️ This pipeline currently runs in **modules-only mode** (no
> containers).

------------------------------------------------------------------------

# 📄 Samplesheet Format

CSV format:

    sample,fastq_1,fastq_2

Example:

    ISO1,/path/to/ISO1_R1.fastq.gz,/path/to/ISO1_R2.fastq.gz
    ISO2,/path/to/ISO2_R1.fastq.gz,/path/to/ISO2_R2.fastq.gz

Absolute FASTQ paths are required.

------------------------------------------------------------------------

# ▶️ Usage

## Run with defaults

``` bash
bash bin/run_pipeline.sh
```

## Specify samplesheet

``` bash
bash bin/run_pipeline.sh --samplesheet config/samplesheet.csv
```

## Specify output directory

``` bash
bash bin/run_pipeline.sh --outdir /scratch/my_run
```

## Specify database location

``` bash
bash bin/run_pipeline.sh --db-path /scratch/my_databases/uva_eskape_2026-01-23
```

## Adjust clustering thresholds

``` bash
bash bin/run_pipeline.sh \
  --round1-threshold 200 \
  --round2-threshold 50
```

------------------------------------------------------------------------

# 🔬 Pipeline Workflow

```{=html}
<details>
```
```{=html}
<summary>
```
`<strong>`{=html}Step 1 --- QC`</strong>`{=html}
```{=html}
</summary>
```
-   TrimGalore\
-   Trimmomatic

Output:

    OUTDIR/qc/<sample>/

```{=html}
</details>
```
```{=html}
<details>
```
```{=html}
<summary>
```
`<strong>`{=html}Step 2 --- Mash Speciation`</strong>`{=html}
```{=html}
</summary>
```
-   mash sketch\
-   mash dist\
-   mash screen

Output:

    OUTDIR/mash/<sample>/

Primary species result:

    <sample>.mash-screen.top3_hits.txt

```{=html}
</details>
```
```{=html}
<details>
```
```{=html}
<summary>
```
`<strong>`{=html}Round 1 --- SNV Calling (All
Isolates)`</strong>`{=html}
```{=html}
</summary>
```
Reference selected to minimize Mash triangle distance.

Outputs:

    OUTDIR/variant_calling_round1/

Includes: - core.aln - core.snp-dists.txt - initial_groups.tsv
```{=html}
</details>
```
```{=html}
<details>
```
```{=html}
<summary>
```
`<strong>`{=html}Round 2 --- Within-Group SNV Calling`</strong>`{=html}
```{=html}
</summary>
```
Runs only for groups with ≥ 2 isolates.

Outputs:

    OUTDIR/variant_calling_round2/groups/GroupX/

Includes: - group-specific ref.fa - core.aln - GroupX.final_groups.tsv
```{=html}
</details>
```

------------------------------------------------------------------------

# 📊 Key Outputs

## ✅ Final Deliverable

    OUTDIR/final_output/final_groups_all.tsv

Columns:

-   Isolate\
-   Species\
-   Primary_Group (Round 1)\
-   Secondary_Group (Round 2)

------------------------------------------------------------------------

## 📈 SNP Distance Matrices

Round 1:

    variant_calling_round1/summary/core.snp-dists.txt

Round 2:

    variant_calling_round2/summary/GroupX.core.snp-dists.txt

Copies are also placed in:

    OUTDIR/final_output/

------------------------------------------------------------------------

# 📁 Output Directory Structure

    OUTDIR/
    ├── qc/
    ├── mash/
    ├── variant_calling_round1/
    ├── variant_calling_round2/
    ├── final_output/
    └── slurm_logs/

------------------------------------------------------------------------

# 🏷 Version

**v0.2.0**\
Two-round SNV calling with configurable thresholds and external database
installation.

