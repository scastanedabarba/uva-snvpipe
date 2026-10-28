#!/usr/bin/env bash
#SBATCH --job-name=snippy_core
# Environment is exported explicitly by the launcher

set -euo pipefail

: "${PIPELINE_ROOT:?ERROR: PIPELINE_ROOT is not set}"
SCRIPT_DIR="${PIPELINE_ROOT}/workflow"
source "${SCRIPT_DIR}/../scripts/setup_env.sh"

OUTROOT="$1"          # variant_calling_round1
CLUSTER_SCRIPT="$2"
THRESHOLD="${3:-200}"
METRICS_SCRIPT="${PIPELINE_ROOT}/scripts/corealn_metrics.py"

[[ -f "$CLUSTER_SCRIPT" ]] || { echo "ERROR: missing cluster script: $CLUSTER_SCRIPT" >&2; exit 1; }
[[ -f "$METRICS_SCRIPT" ]] || { echo "ERROR: missing metrics script: $METRICS_SCRIPT" >&2; exit 1; }

TOP_OUTDIR="$(readlink -f "$(dirname "$OUTROOT")")"

snippy_root="${OUTROOT}/snippy"
summary_dir="${OUTROOT}/summary"
mkdir -p "$summary_dir"

ref_fa="${OUTROOT}/findref/ref.fa"
[[ -s "$ref_fa" ]] || { echo "ERROR: missing ref.fa: $ref_fa" >&2; exit 1; }

mapfile -t dirs < <(find "$snippy_root" -maxdepth 1 -mindepth 1 -type d | sort)
[[ "${#dirs[@]}" -gt 0 ]] || { echo "ERROR: no snippy isolate dirs found in $snippy_root" >&2; exit 1; }

echo "[core] n_isolates=${#dirs[@]}"
echo "[core] running snippy-core ..."

cd "$snippy_root"
snippy-core --ref "$ref_fa" "${dirs[@]}" &> "${snippy_root}/snippy_core.log"

echo "[core] running snp-sites ..."
snp-sites -b -c -o "${snippy_root}/phylo.aln" "${snippy_root}/core.full.aln"

echo "[core] calculating core alignment metrics ..."
python3 "$METRICS_SCRIPT" "${snippy_root}/core.full.aln" \
  --output "${summary_dir}/core_alignment_metrics.tsv"

echo "[core] running snp-dists ..."
snp-dists "${snippy_root}/core.aln" > "${summary_dir}/core.snp-dists.txt"

echo "[core] initial clustering (threshold=${THRESHOLD} SNP) + attaching species ..."
python3 "$CLUSTER_SCRIPT" \
  --snp-dists "${summary_dir}/core.snp-dists.txt" \
  --outdir "$TOP_OUTDIR" \
  --threshold "$THRESHOLD" \
  --output "${summary_dir}/initial_groups.tsv"

echo "[core] outputs:"
echo "  ${snippy_root}/core.aln"
echo "  ${snippy_root}/core.full.aln"
echo "  ${snippy_root}/phylo.aln"
echo "  ${snippy_root}/snippy_core.log"
echo "  ${summary_dir}/core.snp-dists.txt"
echo "  ${summary_dir}/core_alignment_metrics.tsv"
echo "  ${summary_dir}/initial_groups.tsv"

