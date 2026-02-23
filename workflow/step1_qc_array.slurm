#!/usr/bin/env bash
#SBATCH --job-name=qc_reads
#SBATCH --export=NONE

set -euo pipefail

SAMPLESHEET="$1"
OUTROOT="$2"
MODE="$3"

task_id="${SLURM_ARRAY_TASK_ID}"

# ---- pull row for this task_id from CSV ----
# Expected columns: sample, fastq_1, fastq_2
read -r sid r1 r2 < <(
  python3 - <<'PY' "$SAMPLESHEET" "$task_id"
import csv, sys
path=sys.argv[1]; idx=int(sys.argv[2])  # 1-based
with open(path, newline='') as f:
    r=csv.DictReader(f)
    rows=[row for row in r if any((v or "").strip() for v in row.values())]
row=rows[idx-1]
# CHANGE THESE KEYS IF YOUR HEADER DIFFERS:
sid=row["sample"].strip()
r1 =row["fastq_1"].strip()
r2 =row["fastq_2"].strip()
print(sid, r1, r2)
PY
)

echo "[qc] task=${task_id} sample=${sid}"
echo "[qc] r1=${r1}"
echo "[qc] r2=${r2}"
echo "[qc] mode=${MODE}"

sample_out="${OUTROOT}/qc/${sid}"
status_dir="${OUTROOT}/status"
mkdir -p "$sample_out" "$status_dir"
cd "$sample_out"

done_flag="${status_dir}/${sid}.qc.done"
if [[ -f "$done_flag" ]]; then
  echo "[qc] already done -> $done_flag"
  exit 0
fi

if [[ "$MODE" == "modules" ]]; then
  module use /project/amr_services/modulefiles/
  module load trimgalore
  module load trimmomatic/0.39
  module load fastqc/0.11.9
  # fastqstats must be available in PATH in your module setup
elif [[ "$MODE" == "apptainer" ]]; then
  module load apptainer || true
else
  echo "ERROR: MODE must be modules or apptainer (got: $MODE)" >&2
  exit 2
fi

run_cmd() {
  if [[ "$MODE" == "apptainer" ]]; then
    # If you want to use your QC container, uncomment and ensure it contains all tools:
    # apptainer exec docker://ghcr.io/scastanedabarba/uvasnvpipe-qc_reads:v0.1.0 "$@"
    "$@"
  else
    "$@"
  fi
}

# ---- Trim Galore (exact flags from your Nextflow module) ----
run_cmd trim_galore --stringency 5 \
  --trim-n \
  --max_n 5 \
  --length 35 \
  --paired \
  -o . \
  "$r1" "$r2" \
  &> "${sid}.trimgalore.log"

tg_r1_fastq="$(find . -maxdepth 1 -name "*val_1.fq.gz" | head -n 1 || true)"
tg_r2_fastq="$(find . -maxdepth 1 -name "*val_2.fq.gz" | head -n 1 || true)"

if [[ -z "${tg_r1_fastq}" || -z "${tg_r2_fastq}" ]]; then
  echo "ERROR: TrimGalore outputs not found for ${sid}" >&2
  echo "PWD=$(pwd)" >&2
  ls -lah >&2
  exit 1
fi

# ---- Trimmomatic PE (exact parameters; run jar explicitly) ----
run_cmd java -jar /project/amr_services/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 \
  "$tg_r1_fastq" "$tg_r2_fastq" \
  "${sid}_R1.trim.fq.gz" "${sid}_R1.trim.unpaired.fq.gz" \
  "${sid}_R2.trim.fq.gz" "${sid}_R2.trim.unpaired.fq.gz" \
  TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70 \
  &> "${sid}.trimmomatic.log"

cat "${sid}_R1.trim.unpaired.fq.gz" "${sid}_R2.trim.unpaired.fq.gz" > "${sid}_SE.trim.fq.gz"
rm -f "${sid}_R1.trim.unpaired.fq.gz" "${sid}_R2.trim.unpaired.fq.gz"

zcat "${sid}_R1.trim.fq.gz" | /project/amr_services/pipelines/amrpipe/scripts/fastqstats > "${sid}_R1.trim.fastqstats.txt"
zcat "${sid}_R2.trim.fq.gz" | /project/amr_services/pipelines/amrpipe/scripts/fastqstats > "${sid}_R2.trim.fastqstats.txt"

# ---- versions.yml (same as your module) ----
tg_ver="$(trim_galore --version 2>&1 | tr '\n' ' ' | sed -E 's/.*version[[:space:]]+([0-9][0-9.]*) .*/\1/' | head -n 1 || true)"
ca_ver="$(cutadapt --version 2>&1 | head -n 1 || true)"
tr_ver="$(java -jar /project/amr_services/tools/Trimmomatic-0.39/trimmomatic-0.39.jar -version 2>&1 | head -n 1 | sed -E 's/[^0-9]*([0-9.]*).*/\1/' || true)"
jv_ver="$(java -version 2>&1 | head -n 1 | sed -E 's/.*version "([^"]+)".*/\1/' || true)"

cat > versions.yml <<END_VERSIONS
"QC_READS":
  trim_galore: "${tg_ver}"
  cutadapt: "${ca_ver}"
  trimmomatic: "${tr_ver}"
  java: "${jv_ver}"
END_VERSIONS

touch "$done_flag"
echo "[qc] done"

