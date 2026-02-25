#!/usr/bin/env bash
#SBATCH --job-name=r2_snippy
#SBATCH --export=NONE

set -euo pipefail

OUTDIR_ROOT="$(readlink -f "$1")"
ROUND2_OUTDIR="$(readlink -f "$2")"

ROUND1_OUTDIR="${OUTDIR_ROOT}/round1"

groups_tsv="${ROUND2_OUTDIR}/manifest/groups.tsv"
members_dir="${ROUND2_OUTDIR}/manifest/group_members"
groups_dir="${ROUND2_OUTDIR}/groups"

[[ -s "$groups_tsv" ]] || { echo "ERROR: missing groups.tsv: $groups_tsv" >&2; exit 1; }

module load miniforge
set +u
source activate /project/amr_services/.conda/amrpipe-snippy-env
set -u
module load java

echo "[r2_step7] reading groups: $groups_tsv"
echo "[r2_step7] qc reads from: ${ROUND1_OUTDIR}/qc"

while read -r grp nmem; do
  [[ "$grp" == "GroupID" ]] && continue

  gdir="${groups_dir}/${grp}"
  ref_fa="${gdir}/findref/ref.fa"
  [[ -s "$ref_fa" ]] || { echo "ERROR: missing group ref.fa for ${grp}: $ref_fa" >&2; exit 1; }

  done_flag="${ROUND2_OUTDIR}/manifest/${grp}.snippy.done"
  if [[ -f "$done_flag" ]]; then
    echo "[r2_step7] ${grp} snippy already done -> $done_flag"
    continue
  fi

  echo "[r2_step7] group=${grp} n=${nmem}"

  while read -r sid; do
    [[ -z "$sid" ]] && continue

    qc_dir="${ROUND1_OUTDIR}/qc/${sid}"
    r1="${qc_dir}/${sid}_R1.trim.fq.gz"
    r2="${qc_dir}/${sid}_R2.trim.fq.gz"
    [[ -s "$r1" ]] || { echo "ERROR: missing QC R1 for ${sid}: $r1" >&2; exit 1; }
    [[ -s "$r2" ]] || { echo "ERROR: missing QC R2 for ${sid}: $r2" >&2; exit 1; }

    snippy_out="${gdir}/snippy/${sid}"
    mkdir -p "$snippy_out"

    if [[ -f "${snippy_out}/snippy.done" ]]; then
      echo "[r2_step7] ${grp}/${sid} already done"
      continue
    fi

    echo "[r2_step7] running snippy: grp=${grp} sid=${sid}"
    snippy \
      --cpus "${SLURM_CPUS_PER_TASK:-1}" \
      --outdir "$snippy_out" \
      --force \
      --ref "$ref_fa" \
      --R1 "$r1" --R2 "$r2" \
      &> "${snippy_out}/snippy.log"

    touch "${snippy_out}/snippy.done"
  done < "${members_dir}/${grp}.samples.txt"

  touch "$done_flag"
done < "$groups_tsv"

echo "[r2_step7] done"
