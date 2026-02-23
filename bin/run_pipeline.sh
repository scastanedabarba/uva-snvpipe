#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
UVA SNV Pipeline

USAGE:
  run_pipeline.sh [options]

OPTIONS:
  --samplesheet FILE   Path to samplesheet CSV
  --outdir DIR         Output directory (overrides config.sh OUTDIR)
  -h, --help           Show this help message

Samplesheet format (CSV):
  sample,fastq_1,fastq_2

Example:
  run_pipeline.sh --samplesheet config/samplesheet.csv --outdir results/test_run
EOF
}

# -----------------------------
# Parse arguments
# -----------------------------
SAMPLESHEET=""
USER_OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet)
      SAMPLESHEET="$2"
      shift 2
      ;;
    --outdir)
      USER_OUTDIR="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1"
      usage
      exit 1
      ;;
  esac
done

# -----------------------------
# Repo paths
# -----------------------------
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="${REPO_DIR}/bin"
WORKFLOW_DIR="${REPO_DIR}/workflow"
SCRIPTS_DIR="${REPO_DIR}/scripts"
CONFIG_DIR="${REPO_DIR}/config"

# -----------------------------
# Load config defaults
# -----------------------------
source "${CONFIG_DIR}/config.sh"

# -----------------------------
# Override OUTDIR if provided
# -----------------------------
if [[ -n "$USER_OUTDIR" ]]; then
  OUTDIR="$(readlink -f "$USER_OUTDIR")"
fi

# Default samplesheet if not provided
if [[ -z "$SAMPLESHEET" ]]; then
  SAMPLESHEET="${CONFIG_DIR}/samplesheet.csv"
fi

[[ -f "$SAMPLESHEET" ]] || { echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2; exit 1; }

LOGDIR="${OUTDIR}/slurm_logs"

# -----------------------------
# Create output structure
# -----------------------------
mkdir -p "$OUTDIR" "$LOGDIR" \
  "$OUTDIR/qc" \
  "$OUTDIR/mash" \
  "$OUTDIR/findref" \
  "$OUTDIR/snippy" \
  "$OUTDIR/summary" \
  "$OUTDIR/status"

# -----------------------------
# Count samples
# -----------------------------
NSAMPLES=$(
  python3 - <<'PY' "$SAMPLESHEET"
import csv, sys
path=sys.argv[1]
with open(path, newline='') as f:
    r=csv.DictReader(f)
    rows=[row for row in r if any((v or "").strip() for v in row.values())]
print(len(rows))
PY
)

[[ "$NSAMPLES" -gt 0 ]] || { echo "ERROR: No samples found in $SAMPLESHEET" >&2; exit 1; }

echo "[pipeline] samplesheet=$SAMPLESHEET"
echo "[pipeline] nsamples=$NSAMPLES"
echo "[pipeline] outdir=$OUTDIR"
echo "[pipeline] logdir=$LOGDIR"
echo "[pipeline] mode=$MODE"

# -----------------------------
# Step 1: QC
# -----------------------------
QC_JOBID=$(
  sbatch --parsable \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c "$CPUS_QC" --mem "$MEM_QC" -t "$TIME_QC" \
    --array=1-"$NSAMPLES" \
    -J "qc_reads" \
    -o "${LOGDIR}/qc_reads.%A_%a.out" -e "${LOGDIR}/qc_reads.%A_%a.err" \
    "${WORKFLOW_DIR}/step1_qc_array.slurm" \
      "$SAMPLESHEET" "$OUTDIR" "$MODE"
)
echo "[submit] QC array job: $QC_JOBID"

# -----------------------------
# Step 2: MASH
# -----------------------------
MASH_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$QC_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c "$CPUS_MASH" --mem "$MEM_MASH" -t "$TIME_MASH" \
    --array=1-"$NSAMPLES" \
    -J "mash" \
    -o "${LOGDIR}/mash.%A_%a.out" -e "${LOGDIR}/mash.%A_%a.err" \
    "${WORKFLOW_DIR}/step2_mash_array.slurm" \
      "$SAMPLESHEET" "$OUTDIR" "$MODE" "$MASH_DB"
)
echo "[submit] MASH array job: $MASH_JOBID"

# -----------------------------
# Step 3: FINDREF
# -----------------------------
FINDREF_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$MASH_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c 8 --mem "32G" -t "02:00:00" \
    -J "findref" \
    -o "${LOGDIR}/findref.%j.out" -e "${LOGDIR}/findref.%j.err" \
    "${WORKFLOW_DIR}/step3_findref.slurm" \
      "$SAMPLESHEET" "$OUTDIR" "${SCRIPTS_DIR}/parse_mash_triangle.v2.py"
)
echo "[submit] FINDREF job: $FINDREF_JOBID"

# -----------------------------
# Step 4: SNIPPY
# -----------------------------
SNIPPY_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$FINDREF_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c "$CPUS_SNIPPY" --mem "$MEM_SNIPPY" -t "$TIME_SNIPPY" \
    --array=1-"$NSAMPLES" \
    -J "snippy" \
    -o "${LOGDIR}/snippy.%A_%a.out" -e "${LOGDIR}/snippy.%A_%a.err" \
    "${WORKFLOW_DIR}/step4_snippy_array.slurm" \
      "$SAMPLESHEET" "$OUTDIR" "$MODE"
)
echo "[submit] SNIPPY array job: $SNIPPY_JOBID"

# -----------------------------
# Step 5: CORE + CLUSTER
# -----------------------------
CORE_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$SNIPPY_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c 16 --mem "64G" -t "12:00:00" \
    -J "snippy_core" \
    -o "${LOGDIR}/snippy_core.%j.out" -e "${LOGDIR}/snippy_core.%j.err" \
    "${WORKFLOW_DIR}/step5_core.slurm" \
      "$OUTDIR" "${SCRIPTS_DIR}/cluster_isolates.py"
)
echo "[submit] CORE job: $CORE_JOBID"

echo "[pipeline] Submitted:"
echo "  QC      : $QC_JOBID"
echo "  MASH    : $MASH_JOBID"
echo "  FINDREF : $FINDREF_JOBID"
echo "  SNIPPY  : $SNIPPY_JOBID"
echo "  CORE    : $CORE_JOBID"
