#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
UVA SNV Pipeline (SLURM orchestrator)

USAGE:
  run_pipeline.sh [samplesheet.csv]

OPTIONS:
  -h, --help        Show this help message and exit

ARGUMENTS:
  samplesheet.csv   Optional path to a CSV with columns:
                    sample, fastq_1, fastq_2
                    Default: config/samplesheet.csv

NOTES:
  - Configuration is loaded from: config/config.sh
  - This version is modules-only (no containers).
  - Outputs and logs are written to paths defined in config/config.sh (OUTDIR, LOGDIR).

EXAMPLE:
  bin/run_pipeline.sh config/samplesheet.csv
EOF
}

# Handle help flags
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="${REPO_DIR}/bin"
WORKFLOW_DIR="${REPO_DIR}/workflow"
SCRIPTS_DIR="${REPO_DIR}/scripts"
CONFIG_DIR="${REPO_DIR}/config"

# Load pipeline configuration
source "${CONFIG_DIR}/config.sh"

# Default samplesheet location
SAMPLESHEET="${1:-${CONFIG_DIR}/samplesheet.csv}"
[[ -f "$SAMPLESHEET" ]] || { echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2; exit 1; }

mkdir -p "$OUTDIR" "$LOGDIR" \
  "$OUTDIR/qc" "$OUTDIR/mash" "$OUTDIR/findref" "$OUTDIR/snippy" "$OUTDIR/summary" "$OUTDIR/status"

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
echo "[pipeline] repo_dir=$REPO_DIR"
echo "[pipeline] workflow_dir=$WORKFLOW_DIR"
echo "[pipeline] outdir=$OUTDIR"
echo "[pipeline] logdir=$LOGDIR"
echo "[pipeline] mode=$MODE"

# ---------- Step 1: QC array ----------
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

# ---------- Step 2: MASH array (after QC) ----------
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

# ---------- Step 3: FINDREF (single job after MASH) ----------
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

# ---------- Step 4: SNIPPY array (after FINDREF) ----------
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

# ---------- Step 5: CORE + SNP-DISTS (single job after SNIPPY) ----------
CORE_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$SNIPPY_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c 16 --mem "64G" -t "12:00:00" \
    -J "snippy_core" \
    -o "${LOGDIR}/snippy_core.%j.out" -e "${LOGDIR}/snippy_core.%j.err" \
    "${WORKFLOW_DIR}/step5_core.slurm" \
      "$OUTDIR"
)
echo "[submit] CORE job: $CORE_JOBID"

echo "[pipeline] Submitted:"
echo "  QC      : $QC_JOBID"
echo "  MASH    : $MASH_JOBID"
echo "  FINDREF : $FINDREF_JOBID"
echo "  SNIPPY  : $SNIPPY_JOBID"
echo "  CORE    : $CORE_JOBID"
