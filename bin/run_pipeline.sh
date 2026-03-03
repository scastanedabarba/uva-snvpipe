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
  --db-path DIR        Path to database root directory
                       (default: <repo>/databases/uva_eskape_2026-01-23)
  --round1-threshold INT  SNP threshold for initial clustering (default: 200)
  --round2-threshold INT  SNP threshold for final clustering within groups (default: 50)
  -h, --help           Show this help message and exit

SAMPLESHEET FORMAT (CSV):
  sample,fastq_1,fastq_2

OUTPUT LAYOUT:
  <OUTDIR>/qc/
  <OUTDIR>/mash/
  <OUTDIR>/variant_calling_round1/
  <OUTDIR>/variant_calling_round2/
  <OUTDIR>/final_output/
  <OUTDIR>/slurm_logs/
EOF
}

SAMPLESHEET=""
USER_OUTDIR=""
ROUND1_THRESHOLD="200"
ROUND2_THRESHOLD="50"
DB_PATH=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --outdir)      USER_OUTDIR="$2"; shift 2 ;;
    --db-path)     DB_PATH="$2"; shift 2 ;;
    --round1-threshold) ROUND1_THRESHOLD="$2"; shift 2 ;;
    --round2-threshold) ROUND2_THRESHOLD="$2"; shift 2 ;;
    -h|--help)     usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKFLOW_DIR="${REPO_DIR}/workflow"
SCRIPTS_DIR="${REPO_DIR}/scripts"
CONFIG_DIR="${REPO_DIR}/config"

source "${CONFIG_DIR}/config.sh"

DEFAULT_DB_PATH="${REPO_DIR}/databases/uva_eskape_2026-01-23"
if [[ -z "${DB_PATH}" ]]; then
  DB_PATH="$DEFAULT_DB_PATH"
fi
DB_PATH="$(readlink -f "$DB_PATH")"

[[ -d "$DB_PATH" ]] || { echo "ERROR: db-path not found: $DB_PATH" >&2; exit 1; }
[[ -f "$DB_PATH/mash/pathogen_refseq.chrom.k21s50000.msh" ]] || {
  echo "ERROR: missing mash sketch: $DB_PATH/mash/pathogen_refseq.chrom.k21s50000.msh" >&2; exit 1; }
[[ -f "$DB_PATH/mash/pathogen_refseq.chrom.mapping.txt" ]] || {
  echo "ERROR: missing mash mapping: $DB_PATH/mash/pathogen_refseq.chrom.mapping.txt" >&2; exit 1; }
[[ -d "$DB_PATH/chrom_fna" ]] || {
  echo "ERROR: missing chrom_fna/: $DB_PATH/chrom_fna" >&2; exit 1; }

OUTDIR_ROOT="$OUTDIR"
if [[ -n "$USER_OUTDIR" ]]; then
  OUTDIR_ROOT="$(readlink -f "$USER_OUTDIR")"
fi

QC_OUTDIR="${OUTDIR_ROOT}/qc"
MASH_OUTDIR="${OUTDIR_ROOT}/mash"

VC1_OUTDIR="${OUTDIR_ROOT}/variant_calling_round1"
VC2_OUTDIR="${OUTDIR_ROOT}/variant_calling_round2"

FINAL_OUTDIR="${OUTDIR_ROOT}/final_output"
LOGDIR="${OUTDIR_ROOT}/slurm_logs"

if [[ -z "$SAMPLESHEET" ]]; then
  SAMPLESHEET="${CONFIG_DIR}/samplesheet.csv"
fi
[[ -f "$SAMPLESHEET" ]] || { echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2; exit 1; }

mkdir -p "$OUTDIR_ROOT" "$LOGDIR" "$QC_OUTDIR" "$MASH_OUTDIR" "$VC1_OUTDIR" "$VC2_OUTDIR" "$FINAL_OUTDIR"

# round1 dirs
mkdir -p \
  "$VC1_OUTDIR/findref" \
  "$VC1_OUTDIR/snippy" \
  "$VC1_OUTDIR/summary" \
  "$VC1_OUTDIR/status"

# round2 dirs
mkdir -p \
  "$VC2_OUTDIR/manifest" \
  "$VC2_OUTDIR/groups" \
  "$VC2_OUTDIR/summary"

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
echo "[pipeline] db_path=$DB_PATH"
echo "[pipeline] outdir_root=$OUTDIR_ROOT"
echo "[pipeline] qc_dir=$QC_OUTDIR"
echo "[pipeline] mash_dir=$MASH_OUTDIR"
echo "[pipeline] vc1_dir=$VC1_OUTDIR"
echo "[pipeline] vc2_dir=$VC2_OUTDIR"
echo "[pipeline] final_outdir=$FINAL_OUTDIR"
echo "[pipeline] logdir=$LOGDIR"
echo "[pipeline] mode=$MODE"

# ---------- Step 1: QC (writes to OUTDIR_ROOT/qc) ----------
QC_JOBID=$(
  sbatch --parsable \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c "$CPUS_QC" --mem "$MEM_QC" -t "$TIME_QC" \
    --array=1-"$NSAMPLES" \
    -J "qc_reads" \
    -o "${LOGDIR}/step1_qc.%A_%a.out" -e "${LOGDIR}/step1_qc.%A_%a.err" \
    "${WORKFLOW_DIR}/step1_qc_array.slurm" \
      "$SAMPLESHEET" "$OUTDIR_ROOT" "$MODE"
)
echo "[submit] Step1 QC job: $QC_JOBID"

# ---------- Step 2: MASH (reads OUTDIR_ROOT/qc, writes OUTDIR_ROOT/mash) ----------
MASH_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$QC_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c "$CPUS_MASH" --mem "$MEM_MASH" -t "$TIME_MASH" \
    --array=1-"$NSAMPLES" \
    -J "mash" \
    -o "${LOGDIR}/step2_mash.%A_%a.out" -e "${LOGDIR}/step2_mash.%A_%a.err" \
    "${WORKFLOW_DIR}/step2_mash_array.slurm" \
      "$SAMPLESHEET" "$OUTDIR_ROOT" "$MODE" "$DB_PATH"
)
echo "[submit] Step2 MASH job: $MASH_JOBID"

# ---------- Step 3: FINDREF (round1; reads OUTDIR_ROOT/mash, writes VC1/findref) ----------
FINDREF_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$MASH_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c 8 --mem "32G" -t "02:00:00" \
    -J "findref_r1" \
    -o "${LOGDIR}/step3_findref.%j.out" -e "${LOGDIR}/step3_findref.%j.err" \
    "${WORKFLOW_DIR}/step3_findref.slurm" \
      "$SAMPLESHEET" "$VC1_OUTDIR" "${SCRIPTS_DIR}/parse_mash_triangle.v2.py" "$DB_PATH" "$OUTDIR_ROOT"
)
echo "[submit] Step3 FINDREF job: $FINDREF_JOBID"

# ---------- Step 4: SNIPPY (round1; reads OUTDIR_ROOT/qc, writes VC1/snippy) ----------
SNIPPY_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$FINDREF_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c "$CPUS_SNIPPY" --mem "$MEM_SNIPPY" -t "$TIME_SNIPPY" \
    --array=1-"$NSAMPLES" \
    -J "snippy_r1" \
    -o "${LOGDIR}/step4_snippy.%A_%a.out" -e "${LOGDIR}/step4_snippy.%A_%a.err" \
    "${WORKFLOW_DIR}/step4_snippy_array.slurm" \
      "$SAMPLESHEET" "$VC1_OUTDIR" "$MODE" "$OUTDIR_ROOT"
)
echo "[submit] Step4 SNIPPY job: $SNIPPY_JOBID"

# ---------- Step 5: CORE + INITIAL CLUSTERING (threshold=200; mash species from OUTDIR_ROOT/mash) ----------
CORE_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$SNIPPY_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c 16 --mem "64G" -t "12:00:00" \
    -J "core_r1" \
    -o "${LOGDIR}/step5_core.%j.out" -e "${LOGDIR}/step5_core.%j.err" \
    "${WORKFLOW_DIR}/step5_core.slurm" \
      "$VC1_OUTDIR" "${SCRIPTS_DIR}/cluster_isolates.py" "$ROUND1_THRESHOLD"
)
echo "[submit] Step5 CORE+CLUSTER job: $CORE_JOBID"

# ---------- Step 6: ROUND2 FINDREF (prep manifest + per-group findref) ----------
R2_FINDREF_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$CORE_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c 2 --mem "8G" -t "01:00:00" \
    -J "r2_findref" \
    -o "${LOGDIR}/step6_r2_findref.%j.out" -e "${LOGDIR}/step6_r2_findref.%j.err" \
    "${WORKFLOW_DIR}/step6_round2_findref.slurm" \
      "$SAMPLESHEET" "$OUTDIR_ROOT" "$VC1_OUTDIR" "$VC2_OUTDIR" \
      "${SCRIPTS_DIR}/parse_mash_triangle.v2.py" "$DB_PATH"
)
echo "[submit] Step6 ROUND2 FINDREF job: $R2_FINDREF_JOBID"

# ---------- Step 7: ROUND2 SNIPPY (parallel array over ALL isolates; skips non-round2 members) ----------
R2_SNIPPY_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$R2_FINDREF_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c "$CPUS_SNIPPY" --mem "$MEM_SNIPPY" -t "$TIME_SNIPPY" \
    --array=1-"$NSAMPLES" \
    -J "snippy_r2" \
    -o "${LOGDIR}/step7_r2_snippy.%A_%a.out" -e "${LOGDIR}/step7_r2_snippy.%A_%a.err" \
    "${WORKFLOW_DIR}/step7_round2_snippy.slurm" \
      "$SAMPLESHEET" "$OUTDIR_ROOT" "$VC2_OUTDIR" "$MODE"
)
echo "[submit] Step7 ROUND2 SNIPPY job: $R2_SNIPPY_JOBID"

# ---------- Step 8: ROUND2 CORE + FINAL CLUSTERING + MERGE ----------
R2_CORE_JOBID=$(
  sbatch --parsable \
    --dependency=afterok:"$R2_SNIPPY_JOBID" \
    -A "$ACCOUNT" -p "$PARTITION" \
    -c 4 --mem "16G" -t "02:00:00" \
    -J "core_r2" \
    -o "${LOGDIR}/step8_r2_core.%j.out" -e "${LOGDIR}/step8_r2_core.%j.err" \
    "${WORKFLOW_DIR}/step8_round2_core.slurm" \
      "$OUTDIR_ROOT" "$VC2_OUTDIR" "$FINAL_OUTDIR" "${SCRIPTS_DIR}/cluster_isolates.py" "$ROUND2_THRESHOLD"
)
echo "[submit] Step8 ROUND2 CORE+FINAL job: $R2_CORE_JOBID"

echo "[pipeline] Submitted:"
echo "  Step1 QC                 : $QC_JOBID"
echo "  Step2 MASH               : $MASH_JOBID"
echo "  Step3 FINDREF (r1)       : $FINDREF_JOBID"
echo "  Step4 SNIPPY (r1)        : $SNIPPY_JOBID"
echo "  Step5 CORE+INIT (r1)     : $CORE_JOBID"
echo "  Step6 FINDREF (r2)       : $R2_FINDREF_JOBID"
echo "  Step7 SNIPPY (r2)        : $R2_SNIPPY_JOBID"
echo "  Step8 CORE+MERGE (r2)    : $R2_CORE_JOBID"
echo "  Final output dir         : $FINAL_OUTDIR"
