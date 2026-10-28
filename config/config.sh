#!/usr/bin/env bash
# Site-specific settings. Edit before running.
# Absolute path to a Conda environment shared by login and compute nodes.
CONDA_ENV="/sfs/weka/scratch/sgj4qr/snvpipeline/snvenv"

# Reference database installation directory
DB_DIR="/scratch/sgj4qr/snvpipeline/databases/uva_eskape_2026-01-23"

# SLURM settings: leave ACCOUNT or PARTITION empty to use cluster defaults.
ACCOUNT="amr_services_paid"
PARTITION="standard"

# Resources per step (CPU, memory, walltime)
CPUS_QC=1; MEM_QC="8G"; TIME_QC="01:00:00"
CPUS_MASH=4; MEM_MASH="250G"; TIME_MASH="04:00:00"
CPUS_FINDREF=8; MEM_FINDREF="32G"; TIME_FINDREF="02:00:00"
CPUS_SNIPPY=16; MEM_SNIPPY="64G"; TIME_SNIPPY="1:00:00"
CPUS_CORE=16; MEM_CORE="64G"; TIME_CORE="12:00:00"
CPUS_R2_FINDREF=2; MEM_R2_FINDREF="8G"; TIME_R2_FINDREF="01:00:00"
CPUS_R2_CORE=4; MEM_R2_CORE="16G"; TIME_R2_CORE="02:00:00"
