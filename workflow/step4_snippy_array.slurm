#!/usr/bin/env bash
#SBATCH --job-name=snippy
#SBATCH --export=NONE

set -euo pipefail

SAMPLESHEET="$(readlink -f "$1")"
OUTROOT="$(readlink -f "$2")"
MODE="$3"
QCROOT="$(readlink -f "${4:-$(dirname "$OUTROOT")}")"

task_id="${SLURM_ARRAY_TASK_ID}"

read -r sid r1_raw r2_raw < <(
  python3 - <<'PY' "$SAMPLESHEET" "$task_id"
import csv, sys
path=sys.argv[1]; idx=int(sys.argv[2])
with open(path, newline='') as f:
    r=csv.DictReader(f)
    rows=[row for row in r if any((v or "").strip() for v in row.values())]
row=rows[idx-1]
sid=row["sample"].strip()
r1 =row["fastq_1"].strip()
r2 =row["fastq_2"].strip()
print(sid, r1, r2)
PY
)

echo "[snippy] task=${task_id} sample=${sid}"
echo "[snippy] outroot=${OUTROOT}"
echo "[snippy] qcroot=${QCROOT}"

# modules-only for now
if [[ "$MODE" != "modules" ]]; then
  echo "ERROR: MODE must be 'modules' (got: $MODE)" >&2
  exit 2
fi

ref_fa="${OUTROOT}/findref/ref.fa"
[[ -s "$ref_fa" ]] || { echo "ERROR: missing cluster reference fasta: $ref_fa" >&2; exit 1; }

qc_dir="${QCROOT}/qc/${sid}"
snippy_out="${OUTROOT}/snippy/${sid}"
status_dir="${OUTROOT}/status"
mkdir -p "$snippy_out" "$status_dir"

done_flag="${status_dir}/${sid}.snippy.done"
if [[ -f "$done_flag" ]]; then
  echo "[snippy] already done -> $done_flag"
  exit 0
fi

r1="${qc_dir}/${sid}_R1.trim.fq.gz"
r2="${qc_dir}/${sid}_R2.trim.fq.gz"
[[ -s "$r1" ]] || { echo "ERROR: missing QC R1: $r1" >&2; exit 1; }
[[ -s "$r2" ]] || { echo "ERROR: missing QC R2: $r2" >&2; exit 1; }

# Environment setup (Rivanna)
module purge || true
module load miniforge
module load java

# Activate Snippy env
set +u
source activate /project/amr_services/.conda/amrpipe-snippy-env
set -u

command -v snippy >/dev/null 2>&1 || { echo "ERROR: snippy not found in PATH after env activation" >&2; exit 1; }
command -v java   >/dev/null 2>&1 || { echo "ERROR: java not found in PATH" >&2; exit 1; }

echo "[snippy] ref=${ref_fa}"
echo "[snippy] r1=${r1}"
echo "[snippy] r2=${r2}"
echo "[snippy] out=${snippy_out}"

snippy \
  --cpus "${SLURM_CPUS_PER_TASK:-1}" \
  --outdir "$snippy_out" \
  --force \
  --ref "$ref_fa" \
  --R1 "$r1" --R2 "$r2" \
  &> "${snippy_out}/snippy.log"

touch "$done_flag"
echo "[snippy] done -> $done_flag"
