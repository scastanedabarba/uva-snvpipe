#!/usr/bin/env bash
set -euo pipefail

# Cluster
ACCOUNT="amr_services_paid"
PARTITION="standard"

# Resources per step
CPUS_QC=1
MEM_QC="8G"
TIME_QC="01:00:00"

CPUS_MASH=4
MEM_MASH="250G"
TIME_MASH="04:00:00"

CPUS_SNIPPY=16
MEM_SNIPPY="64G"
TIME_SNIPPY="1:00:00"

# Database
MASH_DB="/scratch/sgj4qr/snv_pipeline/snvdb_20260123/mash/pathogen_refseq.chrom.k21s50000.msh"

# Output layout
OUTDIR="$(pwd)/results"
LOGDIR="$(pwd)/slurm_logs"

# Execution mode: "modules" or "apptainer"
# - modules: use module load trimgalore/trimmomatic/fastqstats/mash/etc
# - apptainer: run docker:// containers via apptainer (more reproducible)
MODE="modules"

