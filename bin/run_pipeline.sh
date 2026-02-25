#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
UVA SNV Pipeline

USAGE:
  run_pipeline.sh [options]

OPTIONS:
  --samplesheet FILE   Path to samplesheet CSV (default: config/samplesheet.csv)
  --outdir DIR         Output root directory (overrides config.sh OUTDIR)
  -h, --help           Show this help message and exit

SAMPLESHEET FORMAT (CSV):
  sample,fastq_1,fastq_2

OUTPUT LAYOUT:
  <OUTDIR>/round1/        initial run
  <OUTDIR>/round2/        per-initial-group reruns
  <OUTDIR>/final_output/  merged final results (easy to find)
EOF
}

# -----------------------------
# Parse arguments
# -----------------------------
SAMPLESHEET=""
USER_OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --outdir)      USER_OUTDIR="$2"; shift 2 ;;
    -h|--help)     usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

# -----------------------------
# Repo paths
# -----------------------------
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKFLOW_DIR="${REPO_DIR}/workflow"
SCRIPTS_DIR="${REPO_DIR}/scripts"
CONFIG_DIR="${REPO_DIR}/config"

# -----------------------------
# Load config defaults
# -----------------------------
source "${CONFIG_DIR}/config.sh"

# -----------------------------
# Determine output directories
# -----------------------------
OUTDIR_ROOT="$OUTDIR"
if [[ -n "$USER_OUTDIR" ]]; then
  OUTDIR_ROOT="$(readlink -f "$USER_OUTDIR")"
fi

ROUND1_OUTDIR="${OUTDIR_ROOT}/round1"
ROUND2_OUTDIR="${OUTDIR_ROOT}/round2"
FINAL_OUTDIR="${OUTDIR_ROOT}/final_output"
LOGDIR="${OUTDIR_ROOT}/slurm_logs"

# Default samplesheet
if [[ -z "$SAMPLESHEET" ]]; then
  SAMPLESHEET="${CONFIG_DIR}/samplesheet.csv"
fi
[[ -f "$SAMPLESHEET" ]] || { echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2; exit 1; }

# -----------------------------
# Create output structure (top level + round1 base)
# -----------------------------
mkdir -p "$OUTDIR_ROOT" "$LOGDIR" "$ROUND1_OUTDIR" "$ROUND2_OUTDIR" "$FINAL_OUTDIR"

mkdir -p \
  "$ROUND1_OUTDIR/qc" "$ROUND1_OUTDIR/mash" "$ROUND1_OUTDIR/findref" \
  "$ROUND1_OUTDIR/snippy" "$ROUND1_OUTDIR/summary" "$ROUND1_OUTDIR/status"

mkdir -p \
  "$ROUND2_OUTDIR/clusters" \
  "$ROUND2_OUTDIR/summary"   # may be unused if final_output is used for merged results

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
echo "[pipeline] outdir_root=$OUTDIR_ROOT"
echo "[pipeline] round1_outdir=$ROUND1_OUTDIR"
echo "[pipeline] round2_outdir=$ROUND2_OUTDIR"
echo "[pipeline] final_outdir=$FINAL_OUTDIR"
echo "[pipeline] logdir=$LOGDIR"
echo "[pipeline] mode=$MODE"

# ---------------------------------------------------------------------
# ROUND 1 (same as before, but OUTROOT=ROUND1_OUTDIR)
# ---------------------------------------------------------------------

# ---------- Step 1: QC ----------
QC_JOBID=$(
  sbatch --parsable \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c "$CPUS_QC" --mem "$MEM_QC" -t "$TIME_QC" \
    --array=1-"$NSAMPLES" \
    -J "qc_reads_r1" \
    -o "${LOGDIR}/step1_qc.%A_%a.out" -e "${LOGDIR}/step1_qc.%A_%a.err" \
    "${WORKFLOW_DIR}/step1_qc_array.slurm" \
      "$SAMPLESHEET" "$ROUND1_OUTDIR" "$MODE"
)
echo "[submit] Step1 QC job: $QC_JOBID"

# ---------- Step 2: MASH ----------
MASH_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$QC_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c "$CPUS_MASH" --mem "$MEM_MASH" -t "$TIME_MASH" \
    --array=1-"$NSAMPLES" \
    -J "mash_r1" \
    -o "${LOGDIR}/step2_mash.%A_%a.out" -e "${LOGDIR}/step2_mash.%A_%a.err" \
    "${WORKFLOW_DIR}/step2_mash_array.slurm" \
      "$SAMPLESHEET" "$ROUND1_OUTDIR" "$MODE" "$MASH_DB"
)
echo "[submit] Step2 MASH job: $MASH_JOBID"

# ---------- Step 3: FINDREF ----------
FINDREF_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$MASH_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c 8 --mem "32G" -t "02:00:00" \
    -J "findref_r1" \
    -o "${LOGDIR}/step3_findref.%j.out" -e "${LOGDIR}/step3_findref.%j.err" \
    "${WORKFLOW_DIR}/step3_findref.slurm" \
      "$SAMPLESHEET" "$ROUND1_OUTDIR" "${SCRIPTS_DIR}/parse_mash_triangle.v2.py"
)
echo "[submit] Step3 FINDREF job: $FINDREF_JOBID"

# ---------- Step 4: SNIPPY ----------
SNIPPY_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$FINDREF_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c "$CPUS_SNIPPY" --mem "$MEM_SNIPPY" -t "$TIME_SNIPPY" \
    --array=1-"$NSAMPLES" \
    -J "snippy_r1" \
    -o "${LOGDIR}/step4_snippy.%A_%a.out" -e "${LOGDIR}/step4_snippy.%A_%a.err" \
    "${WORKFLOW_DIR}/step4_snippy_array.slurm" \
      "$SAMPLESHEET" "$ROUND1_OUTDIR" "$MODE"
)
echo "[submit] Step4 SNIPPY job: $SNIPPY_JOBID"

# ---------- Step 5: CORE + INITIAL CLUSTERING ----------
# This will be updated to:
#   - produce initial_groups.tsv with threshold=200
#   - write core artifacts into ROUND1_OUTDIR/snippy or ROUND1_OUTDIR/core depending on your earlier structure
CORE_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$SNIPPY_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c 16 --mem "64G" -t "12:00:00" \
    -J "core_r1" \
    -o "${LOGDIR}/step5_core.%j.out" -e "${LOGDIR}/step5_core.%j.err" \
    "${WORKFLOW_DIR}/step5_core.slurm" \
      "$ROUND1_OUTDIR" "${SCRIPTS_DIR}/cluster_isolates.py"
)
echo "[submit] Step5 CORE+CLUSTER job: $CORE_JOBID"

# ---------------------------------------------------------------------
# ROUND 2 (sequential steps, like round1)
# Step6 will prep round2 dir structure AND submit per-group findref jobs.
# Step7 submits per-group snippy jobs.
# Step8 submits per-group core+final clustering jobs AND merges to final_output.
# ---------------------------------------------------------------------

# ---------- Step 6: ROUND2 FINDREF (prep + submit per-group findref) ----------
R2_FINDREF_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$CORE_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c 2 --mem "8G" -t "01:00:00" \
    -J "r2_findref" \
    -o "${LOGDIR}/step6_r2_findref.%j.out" -e "${LOGDIR}/step6_r2_findref.%j.err" \
    "${WORKFLOW_DIR}/step6_round2_findref.slurm" \
      "$SAMPLESHEET" "$OUTDIR_ROOT" "$ROUND1_OUTDIR" "$ROUND2_OUTDIR" \
      "${SCRIPTS_DIR}/parse_mash_triangle.v2.py"
)
echo "[submit] Step6 ROUND2 FINDREF job: $R2_FINDREF_JOBID"

# ---------- Step 7: ROUND2 SNIPPY (submit per-group snippy arrays) ----------
R2_SNIPPY_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$R2_FINDREF_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c 2 --mem "8G" -t "01:00:00" \
    -J "r2_snippy" \
    -o "${LOGDIR}/step7_r2_snippy.%j.out" -e "${LOGDIR}/step7_r2_snippy.%j.err" \
    "${WORKFLOW_DIR}/step7_round2_snippy.slurm" \
      "$OUTDIR_ROOT" "$ROUND2_OUTDIR"
)
echo "[submit] Step7 ROUND2 SNIPPY job: $R2_SNIPPY_JOBID"

# ---------- Step 8: ROUND2 CORE + FINAL CLUSTERING (submit per-group core jobs + merge final_output) ----------
R2_CORE_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$R2_SNIPPY_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c 2 --mem "8G" -t "01:00:00" \
    -J "r2_core" \
    -o "${LOGDIR}/step8_r2_core.%j.out" -e "${LOGDIR}/step8_r2_core.%j.err" \
    "${WORKFLOW_DIR}/step8_round2_core.slurm" \
      "$OUTDIR_ROOT" "$ROUND2_OUTDIR" "$FINAL_OUTDIR" "${SCRIPTS_DIR}/cluster_isolates.py"
)
echo "[submit] Step8 ROUND2 CORE+FINAL job: $R2_CORE_JOBID"

echo "[pipeline] Submitted:"
echo "  Step1 QC                 : $QC_JOBID"
echo "  Step2 MASH               : $MASH_JOBID"
echo "  Step3 FINDREF            : $FINDREF_JOBID"
echo "  Step4 SNIPPY             : $SNIPPY_JOBID"
echo "  Step5 CORE+INIT CLUSTER  : $CORE_JOBID"
echo "  Step6 R2 FINDREF         : $R2_FINDREF_JOBID"
echo "  Step7 R2 SNIPPY          : $R2_SNIPPY_JOBID"
echo "  Step8 R2 CORE+FINAL      : $R2_CORE_JOBID"
echo "  Final output dir         : $FINAL_OUTDIR"
