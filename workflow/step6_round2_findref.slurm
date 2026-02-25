#!/usr/bin/env bash
#SBATCH --job-name=r2_findref
#SBATCH --export=NONE

set -euo pipefail

SAMPLESHEET="$(readlink -f "$1")"
OUTDIR_ROOT="$(readlink -f "$2")"
ROUND1_OUTDIR="$(readlink -f "$3")"
ROUND2_OUTDIR="$(readlink -f "$4")"
PARSE_TRI="$(readlink -f "$5")"

# ---- Required round1 inputs ----
INITIAL_GROUPS="${ROUND1_OUTDIR}/summary/initial_groups.tsv"
MAP_TSV="/scratch/sgj4qr/snv_pipeline/snvdb_20260123/mash/pathogen_refseq.chrom.mapping.txt"
CHROM_FNA_DIR="/scratch/sgj4qr/snv_pipeline/snvdb_20260123/chrom_fna/"

[[ -f "$SAMPLESHEET" ]]      || { echo "ERROR: missing samplesheet: $SAMPLESHEET" >&2; exit 1; }
[[ -s "$INITIAL_GROUPS" ]]   || { echo "ERROR: missing initial_groups.tsv: $INITIAL_GROUPS" >&2; exit 1; }
[[ -f "$PARSE_TRI" ]]        || { echo "ERROR: missing parse script: $PARSE_TRI" >&2; exit 1; }
[[ -s "$MAP_TSV" ]]          || { echo "ERROR: missing mapping file: $MAP_TSV" >&2; exit 1; }
[[ -d "$CHROM_FNA_DIR" ]]    || { echo "ERROR: missing chrom_fna dir: $CHROM_FNA_DIR" >&2; exit 1; }

module use /project/amr_services/modulefiles/
module load mash/2.3

manifest_dir="${ROUND2_OUTDIR}/manifest"
members_dir="${manifest_dir}/group_members"
groups_dir="${ROUND2_OUTDIR}/groups"

mkdir -p "$manifest_dir" "$members_dir" "$groups_dir"

groups_tsv="${manifest_dir}/groups.tsv"
singletons_tsv="${manifest_dir}/singletons.tsv"

echo -e "GroupID\tn_members" > "$groups_tsv"
echo -e "Isolate\tRound1Group" > "$singletons_tsv"

echo "[r2_step6] reading round1 initial groups: $INITIAL_GROUPS"
echo "[r2_step6] writing manifest: $manifest_dir"

# ---- Build group membership files (exclude Unrelated and size==1) ----
python3 - <<'PY' "$INITIAL_GROUPS" "$groups_tsv" "$members_dir" "$singletons_tsv"
import csv, sys
from collections import defaultdict

initial_groups, groups_tsv, members_dir, singletons_tsv = sys.argv[1:5]

# Read initial_groups.tsv: Isolate Species Group
groups = defaultdict(list)
rows = []
with open(initial_groups, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        iso = row["Isolate"].strip()
        grp = row["Group"].strip()
        rows.append((iso, grp))
        groups[grp].append(iso)

# Determine which groups to keep for round2
kept = []
for grp, members in groups.items():
    if grp == "Unrelated":
        continue
    if len(members) < 2:
        continue
    kept.append((grp, len(members)))

kept.sort(key=lambda x: (int(x[0].replace("Group","")) if x[0].startswith("Group") else 10**9, x[0]))

# Write groups.tsv
with open(groups_tsv, "a") as out:
    for grp, n in kept:
        out.write(f"{grp}\t{n}\n")

# Write group member files
for grp, _n in kept:
    path = f"{members_dir}/{grp}.samples.txt"
    with open(path, "w") as out:
        for iso in sorted(groups[grp]):
            out.write(iso + "\n")

# Write singletons: Unrelated + any size==1 groups
with open(singletons_tsv, "a") as out:
    for iso, grp in sorted(rows):
        if grp == "Unrelated" or len(groups[grp]) < 2:
            out.write(f"{iso}\t{grp}\n")
PY

# ---- Create per-group directory skeleton + per-group samplesheet.csv ----
echo "[r2_step6] creating per-group directories under: $groups_dir"
while read -r grp nmem; do
  [[ "$grp" == "GroupID" ]] && continue
  gdir="${groups_dir}/${grp}"
  mkdir -p \
    "${gdir}/findref" \
    "${gdir}/snippy" \
    "${gdir}/core" \
    "${gdir}/summary"

  # build group samplesheet.csv from global samplesheet
  python3 - <<'PY' "$SAMPLESHEET" "${members_dir}/${grp}.samples.txt" "${gdir}/samplesheet.csv"
import csv, sys
samplesheet, members_file, out_csv = sys.argv[1:4]
members=set(x.strip() for x in open(members_file) if x.strip())

with open(samplesheet, newline="") as f:
    r=csv.DictReader(f)
    rows=[row for row in r if row.get("sample","").strip() in members]

if not rows:
    raise SystemExit(f"ERROR: no rows found for group members in {members_file}")

with open(out_csv, "w", newline="") as out:
    w=csv.DictWriter(out, fieldnames=r.fieldnames)
    w.writeheader()
    for row in rows:
        w.writerow(row)
PY
done < "$groups_tsv"

# ---- Run findref within each group (using round1 mash outputs) ----
echo "[r2_step6] running findref for each round2 group ..."
while read -r grp nmem; do
  [[ "$grp" == "GroupID" ]] && continue

  gdir="${groups_dir}/${grp}"
  gfind="${gdir}/findref"
  gmsh="${gfind}/msh"
  mkdir -p "$gmsh"

  done_flag="${ROUND2_OUTDIR}/manifest/${grp}.findref.done"
  if [[ -f "$done_flag" ]]; then
    echo "[r2_step6] ${grp} findref already done -> $done_flag"
    continue
  fi

  echo "[r2_step6] group=${grp} n=${nmem}"

  # sample IDs
  mapfile -t samples < "${members_dir}/${grp}.samples.txt"

  # 1) copy isolate sketches from round1 mash
  for sid in "${samples[@]}"; do
    iso_msh="${ROUND1_OUTDIR}/mash/${sid}/${sid}.reads.fq.gz.msh"
    [[ -s "$iso_msh" ]] || { echo "ERROR: missing isolate sketch for ${sid}: $iso_msh" >&2; exit 1; }
    cp -f "$iso_msh" "$gmsh/"
  done

  # 2) sketch each group member's top-hit reference (from round1 mash-ref) into gmsh/
  for sid in "${samples[@]}"; do
    mash_ref="${ROUND1_OUTDIR}/mash/${sid}/${sid}.mash-ref.txt"
    [[ -s "$mash_ref" ]] || { echo "ERROR: missing mash-ref for ${sid}: $mash_ref" >&2; exit 1; }

    ref_seqid="$(awk 'NR==1{print $2}' "$mash_ref")"
    [[ -n "$ref_seqid" ]] || { echo "ERROR: empty ref_seqid for $sid in $mash_ref" >&2; exit 1; }

    ref_accid="${ref_seqid%%|*}"
    reffa="$(find "$CHROM_FNA_DIR" -maxdepth 1 -name "${ref_accid}*.fna" | head -n 1 || true)"
    if [[ -z "$reffa" ]]; then
      echo "[r2_step6] WARN: no chrom_fna fasta found for ${ref_accid}; skipping"
      continue
    fi

    cp -f "$reffa" "$gmsh/"
    local_fa="${gmsh}/$(basename "$reffa")"
    mash sketch -m 2 -k 21 -s 50000 "$local_fa"
    rm -f "$local_fa"
  done

  # 3) mash triangle within group
  find "$gmsh" -name "*.msh" | sort > "${gmsh}/msh_linelist"
  mash triangle -p "${SLURM_CPUS_PER_TASK:-2}" -l "${gmsh}/msh_linelist" > "${gfind}/mash-triangle.output.txt"

  # 4) parse to closest ref label
  python3 "$PARSE_TRI" -t "${gfind}/mash-triangle.output.txt" -o "${gfind}/mash-closest-ref.txt"

  ref_label="$(awk 'NR==1{print $1; exit}' "${gfind}/mash-closest-ref.txt")"
  [[ -n "$ref_label" ]] || { echo "ERROR: empty mash-closest-ref.txt for ${grp}" >&2; exit 1; }

  # 5) locate fasta and write ref.fa
  if [[ "$ref_label" == *.fna ]]; then
    ref_fa="${CHROM_FNA_DIR%/}/${ref_label}"
  else
    ref_fa="$(find "$CHROM_FNA_DIR" -maxdepth 1 -name "${ref_label}*.fna" | head -n 1 || true)"
  fi

  [[ -s "$ref_fa" ]] || { echo "ERROR: could not locate fasta for ref_label=${ref_label} (grp=${grp})" >&2; exit 1; }

  cp -f "$ref_fa" "${gfind}/ref.fa"
  echo "$ref_label" > "${gfind}/ref_accid.txt"

  touch "$done_flag"
  echo "[r2_step6] ${grp} chose ref=${ref_label}"
done < "$groups_tsv"

echo "[r2_step6] done"
