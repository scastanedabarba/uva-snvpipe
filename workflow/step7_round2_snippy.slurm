#!/usr/bin/env bash
#SBATCH --job-name=r2_snippy
#SBATCH --export=NONE

set -euo pipefail

SAMPLESHEET="$1"
OUTDIR_ROOT="$(readlink -f "$2")"
ROUND2_OUTDIR="$(readlink -f "$3")"
MODE="$4"

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

# modules-only
if [[ "$MODE" != "modules" ]]; then
  echo "ERROR: MODE must be 'modules' (got: $MODE)" >&2
  exit 2
fi

iso2grp="${ROUND2_OUTDIR}/manifest/isolate_to_group.tsv"
[[ -s "$iso2grp" ]] || { echo "ERROR: missing isolate_to_group.tsv: $iso2grp" >&2; exit 1; }

# lookup group for this isolate; empty => not in round2
grp="$(awk -v s="$sid" 'BEGIN{FS="\t"} NR>1 && $1==s {print $2; exit}' "$iso2grp")"
if [[ -z "${grp}" ]]; then
  echo "[r2_step7] sample=${sid} not in any round2 group -> skipping"
  exit 0
fi

gdir="${ROUND2_OUTDIR}/groups/${grp}"
ref_fa="${gdir}/findref/ref.fa"
[[ -s "$ref_fa" ]] || { echo "ERROR: missing group ref.fa for ${grp}: $ref_fa" >&2; exit 1; }

qc_dir="${OUTDIR_ROOT}/qc/${sid}"
r1="${qc_dir}/${sid}_R1.trim.fq.gz"
r2="${qc_dir}/${sid}_R2.trim.fq.gz"
[[ -s "$r1" ]] || { echo "ERROR: missing QC R1 for ${sid}: $r1" >&2; exit 1; }
[[ -s "$r2" ]] || { echo "ERROR: missing QC R2 for ${sid}: $r2" >&2; exit 1; }

snippy_out="${gdir}/snippy/${sid}"
mkdir -p "$snippy_out"

if [[ -f "${snippy_out}/snippy.done" ]]; then
  echo "[r2_step7] ${grp}/${sid} already done"
  exit 0
fi

module purge || true
module load miniforge
module load java
set +u
source activate /project/amr_services/.conda/amrpipe-snippy-env
set -u

echo "[r2_step7] grp=${grp} sid=${sid}"
echo "[r2_step7] ref=${ref_fa}"
echo "[r2_step7] r1=${r1}"
echo "[r2_step7] r2=${r2}"
echo "[r2_step7] out=${snippy_out}"

snippy \
  --cpus "${SLURM_CPUS_PER_TASK:-1}" \
  --outdir "$snippy_out" \
  --force \
  --ref "$ref_fa" \
  --R1 "$r1" --R2 "$r2" \
  &> "${snippy_out}/snippy.log"

touch "${snippy_out}/snippy.done"
echo "[r2_step7] done"
