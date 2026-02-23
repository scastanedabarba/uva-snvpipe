#!/usr/bin/env bash
#SBATCH --job-name=mash
#SBATCH --export=NONE

set -euo pipefail

SAMPLESHEET="$1"
OUTROOT="$2"
MODE="$3"
MASHREF="$4"          # same as MASH_DB in config, but keeping your original var name

task_id="${SLURM_ARRAY_TASK_ID}"

read -r ACCID R1FASTQ R2FASTQ < <(
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

THREADS="${SLURM_CPUS_PER_TASK:-1}"

echo "[mash] task=${task_id} sample=${ACCID}"
echo "[mash] mode=${MODE}"
echo "[mash] mash_ref=${MASHREF}"
echo "[mash] threads=${THREADS}"

qc_dir="${OUTROOT}/qc/${ACCID}"
OUTDIR="${OUTROOT}/mash/${ACCID}"
status_dir="${OUTROOT}/status"
mkdir -p "$OUTDIR" "$status_dir"
cd "$OUTDIR"

done_flag="${status_dir}/${ACCID}.mash.done"
if [[ -f "$done_flag" ]]; then
  echo "[mash] already done -> $done_flag"
  exit 0
fi

# modules-only for now
if [[ "$MODE" != "modules" ]]; then
  echo "ERROR: MODE must be 'modules' (got: $MODE)" >&2
  exit 2
fi

module use /project/amr_services/modulefiles/
module load mash/2.3

# -----------------------------
# Locate QC reads robustly
# Prefer final Trimmomatic outputs, else fall back to TrimGalore val_1/val_2
# -----------------------------
R1FASTQ="${qc_dir}/${ACCID}_R1.trim.fq.gz"
R2FASTQ="${qc_dir}/${ACCID}_R2.trim.fq.gz"

if [[ ! -s "$R1FASTQ" || ! -s "$R2FASTQ" ]]; then
  tg_r1="$(find "$qc_dir" -maxdepth 1 -name "*val_1.fq.gz" | head -n 1 || true)"
  tg_r2="$(find "$qc_dir" -maxdepth 1 -name "*val_2.fq.gz" | head -n 1 || true)"
  if [[ -n "${tg_r1}" && -n "${tg_r2}" && -s "${tg_r1}" && -s "${tg_r2}" ]]; then
    echo "[mash] WARN: Trimmomatic outputs not found; using TrimGalore val outputs instead"
    R1FASTQ="$tg_r1"
    R2FASTQ="$tg_r2"
  else
    echo "ERROR: Could not find QC reads for ${ACCID}" >&2
    echo "Tried:" >&2
    echo "  $R1FASTQ" >&2
    echo "  $R2FASTQ" >&2
    echo "Also searched for *val_1.fq.gz/*val_2.fq.gz in: $qc_dir" >&2
    ls -lah "$qc_dir" >&2 || true
    exit 1
  fi
fi

echo "[mash] R1FASTQ=${R1FASTQ}"
echo "[mash] R2FASTQ=${R2FASTQ}"

# mapping file sits next to the .msh in your DB directory
MASHREFMAP="$(dirname "$MASHREF")/pathogen_refseq.chrom.mapping.txt"
if [[ ! -s "$MASHREFMAP" ]]; then
  echo "ERROR: Missing MASHREFMAP: $MASHREFMAP" >&2
  exit 1
fi

echo -e "Determine closest species match using MASH distance ... \n"

# EXACT original structure + filenames (gzip-safe)
zcat "$R1FASTQ" "$R2FASTQ" | gzip -c > "${OUTDIR}/${ACCID}.reads.fq.gz"
mash sketch -r -m 2 -k 21 -s 50000 "${OUTDIR}/${ACCID}.reads.fq.gz"

# mash dist
mash dist -p "${THREADS}" "${OUTDIR}/${ACCID}.reads.fq.gz.msh" "${MASHREF}" > "${OUTDIR}/${ACCID}.mash-dist.output.txt"

### get top hit from mash-dist  (guard SIGPIPE from sort|head under pipefail)
set +o pipefail
refaccid=$(sort -gk3 "${OUTDIR}/${ACCID}.mash-dist.output.txt" | head -1 | cut -f 2)
set -o pipefail

refstr=$(awk -v var="$refaccid" '$1 == var' "$MASHREFMAP" | cut -f 2)
echo -e "${ACCID}\t${refaccid}\t${refstr}" > "${OUTDIR}/${ACCID}.mash-ref.txt"

# mash screen
mash screen -p "${THREADS}" -w "${MASHREF}" "${OUTDIR}/${ACCID}.reads.fq.gz" | sort -gr > "${OUTDIR}/${ACCID}.mash-screen.output.txt"

### get top3 hits from screening
head -3 "${OUTDIR}/${ACCID}.mash-screen.output.txt" > "${OUTDIR}/${ACCID}.mash-screen.top3_hits.txt"

touch "$done_flag"
echo "[mash] done -> $done_flag"
