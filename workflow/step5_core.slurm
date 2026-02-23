#!/usr/bin/env bash
#SBATCH --job-name=snippy_core
#SBATCH --export=NONE

set -euo pipefail

OUTROOT="$1"

module load miniforge
set +u
source activate /project/amr_services/.conda/amrpipe-snippy-env
set -u
module use /project/amr_services/modulefiles/
module load snp-dists || true

snippy_root="${OUTROOT}/snippy"
summary_dir="${OUTROOT}/summary"
mkdir -p "$summary_dir"

ref_fa="${OUTROOT}/findref/ref.fa"
[[ -s "$ref_fa" ]] || { echo "ERROR: missing ref.fa: $ref_fa" >&2; exit 1; }

# All per-isolate snippy dirs
mapfile -t dirs < <(find "$snippy_root" -maxdepth 1 -mindepth 1 -type d | sort)
[[ "${#dirs[@]}" -gt 0 ]] || { echo "ERROR: no snippy isolate dirs found in $snippy_root" >&2; exit 1; }

echo "[core] n_isolates=${#dirs[@]}"
echo "[core] running snippy-core ..."

# IMPORTANT: run snippy-core from within snippy_root so core.* outputs land in OUTROOT/snippy/
cd "$snippy_root"
snippy-core --ref "$ref_fa" "${dirs[@]}" &> "${summary_dir}/snippy_core.log"

echo "[core] running snp-sites ..."
snp-sites -b -c -o "${summary_dir}/phylo.aln" "${snippy_root}/core.full.aln"

echo "[core] running snp-dists ..."
snp-dists "${snippy_root}/core.aln" > "${summary_dir}/core.snp-dists.txt"

echo "[core] outputs:"
echo "  ${snippy_root}/core.aln"
echo "  ${snippy_root}/core.full.aln"
echo "  ${summary_dir}/phylo.aln"
echo "  ${summary_dir}/core.snp-dists.txt"
echo "  ${summary_dir}/snippy_core.log"
