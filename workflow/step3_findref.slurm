#!/usr/bin/env bash
#SBATCH --job-name=findref
#SBATCH --export=NONE

set -euo pipefail

SAMPLESHEET="$(readlink -f "$1")"
OUTROOT="$2"
PARSE_TRI="$3"

# NEW: MASHROOT is the directory that contains mash/<isolate>/... (default: parent of OUTROOT)
# Backward-compatible heuristic:
#   - if $4 looks like a mapping file, treat it as MAP_TSV (old calling style)
#   - otherwise treat $4 as MASHROOT
ARG4="${4:-}"
if [[ -n "$ARG4" && "$(basename "$ARG4")" == *mapping*.txt ]]; then
  MASHROOT="$(readlink -f "$(dirname "$OUTROOT")")"
  MAP_TSV="$(readlink -f "$ARG4")"
  CHROM_FNA_DIR="${5:-/scratch/sgj4qr/snv_pipeline/snvdb_20260123/chrom_fna/}"
else
  MASHROOT="$(readlink -f "${ARG4:-$(dirname "$OUTROOT")}")"
  MAP_TSV="${5:-/scratch/sgj4qr/snv_pipeline/snvdb_20260123/mash/pathogen_refseq.chrom.mapping.txt}"
  CHROM_FNA_DIR="${6:-/scratch/sgj4qr/snv_pipeline/snvdb_20260123/chrom_fna/}"
fi

[[ -f "$MAP_TSV" ]]       || { echo "ERROR: missing mapping file: $MAP_TSV" >&2; exit 1; }
[[ -d "$CHROM_FNA_DIR" ]] || { echo "ERROR: missing chrom_fna dir: $CHROM_FNA_DIR" >&2; exit 1; }
[[ -f "$PARSE_TRI" ]]     || { echo "ERROR: missing parse script: $PARSE_TRI" >&2; exit 1; }

module use /project/amr_services/modulefiles/
module load mash/2.3

findref_dir="${OUTROOT}/findref"
msh_dir="${findref_dir}/msh"
mkdir -p "$msh_dir"
cd "$findref_dir"

# Build a simple list of sample IDs from samplesheet.csv
samples=($(python3 - <<'PY' "$SAMPLESHEET"
import csv, sys
path=sys.argv[1]
with open(path, newline='') as f:
    r=csv.DictReader(f)
    rows=[row for row in r if any((v or "").strip() for v in row.values())]
print("\n".join([row["sample"].strip() for row in rows]))
PY
))

echo "[findref] n_samples=${#samples[@]}"
echo "[findref] mashroot=${MASHROOT}"
echo "[findref] outroot=${OUTROOT}"
echo "[findref] copying isolate sketches + sketching candidate refs"

# 1) Copy isolate sketches into msh/
for sid in "${samples[@]}"; do
  iso_msh="${MASHROOT}/mash/${sid}/${sid}.reads.fq.gz.msh"
  [[ -s "$iso_msh" ]] || { echo "ERROR: missing isolate sketch: $iso_msh" >&2; exit 1; }
  cp -f "$iso_msh" "$msh_dir/"
done

# 2) For each isolate: read its top-hit ref_seqid (from mash-ref) and sketch that reference
for sid in "${samples[@]}"; do
  mash_ref="${MASHROOT}/mash/${sid}/${sid}.mash-ref.txt"
  [[ -s "$mash_ref" ]] || { echo "ERROR: missing mash-ref: $mash_ref" >&2; exit 1; }

  ref_seqid="$(awk 'NR==1{print $2}' "$mash_ref")"
  [[ -n "$ref_seqid" ]] || { echo "ERROR: empty ref_seqid for $sid in $mash_ref" >&2; exit 1; }

  ref_accid="${ref_seqid%%|*}"
  [[ -n "$ref_accid" ]] || { echo "ERROR: could not parse ref_accid from ref_seqid=$ref_seqid (sample=$sid)" >&2; exit 1; }

  reffa="$(find "$CHROM_FNA_DIR" -name "${ref_accid}*.fna" | head -n 1 || true)"
  if [[ -z "${reffa}" ]]; then
    echo "WARN: No chrom_fna fasta found for ${ref_accid} (from ${ref_seqid}); skipping"
    continue
  fi

  cp -f "$reffa" "$msh_dir/"
  local_fa="${msh_dir}/$(basename "$reffa")"

  mash sketch -m 2 -k 21 -s 50000 "$local_fa"
  rm -f "$local_fa"
done

# 3) Triangle across all msh
find "$msh_dir" -name "*.msh" | sort > "${msh_dir}/msh_linelist"
mash triangle -p "${SLURM_CPUS_PER_TASK:-8}" -l "${msh_dir}/msh_linelist" > "${findref_dir}/mash-triangle.output.txt"

# 4) Parse to closest-ref
python3 "$PARSE_TRI" -t "${findref_dir}/mash-triangle.output.txt" -o "${findref_dir}/mash-closest-ref.txt"
echo "[findref] wrote: ${findref_dir}/mash-closest-ref.txt"

# 5) Extract chosen reference label and prepare ref.fa for snippy
ref_label="$(awk 'NR==1{print $1; exit}' "${findref_dir}/mash-closest-ref.txt")"
[[ -n "$ref_label" ]] || { echo "ERROR: empty mash-closest-ref.txt" >&2; exit 1; }

echo "$ref_label" > "${findref_dir}/ref_accid.txt"

if [[ "$ref_label" == *.fna ]]; then
  ref_fa="${CHROM_FNA_DIR%/}/${ref_label}"
else
  ref_fa="$(find "$CHROM_FNA_DIR" -maxdepth 1 -name "${ref_label}*.fna" | head -n 1 || true)"
fi

[[ -s "$ref_fa" ]] || {
  echo "ERROR: could not locate fasta for chosen ref_label=$ref_label" >&2
  echo "Tried ref_fa=$ref_fa" >&2
  echo "CHROM_FNA_DIR=$CHROM_FNA_DIR" >&2
  exit 1
}

cp -f "$ref_fa" "${findref_dir}/ref.fa"
echo "[findref] chosen ref: $ref_label"
echo "[findref] ref_fa: $ref_fa"
echo "[findref] wrote: ${findref_dir}/ref.fa"
