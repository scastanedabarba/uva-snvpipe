#!/usr/bin/env bash
#SBATCH --job-name=r2_core
#SBATCH --export=NONE

set -euo pipefail

OUTDIR_ROOT="$(readlink -f "$1")"
ROUND2_OUTDIR="$(readlink -f "$2")"
FINAL_OUTDIR="$(readlink -f "$3")"
CLUSTER_SCRIPT="$(readlink -f "$4")"

ROUND1_OUTDIR="${OUTDIR_ROOT}/round1"

groups_tsv="${ROUND2_OUTDIR}/manifest/groups.tsv"
members_dir="${ROUND2_OUTDIR}/manifest/group_members"
singletons_tsv="${ROUND2_OUTDIR}/manifest/singletons.tsv"
groups_dir="${ROUND2_OUTDIR}/groups"

[[ -s "$groups_tsv" ]] || { echo "ERROR: missing groups.tsv: $groups_tsv" >&2; exit 1; }
[[ -f "$CLUSTER_SCRIPT" ]] || { echo "ERROR: missing cluster script: $CLUSTER_SCRIPT" >&2; exit 1; }

module load miniforge
set +u
source activate /project/amr_services/.conda/amrpipe-snippy-env
set -u
module use /project/amr_services/modulefiles/
module load snp-dists || true

mkdir -p "$FINAL_OUTDIR"

echo "[r2_step8] running per-group core + final clustering ..."
while read -r grp nmem; do
  [[ "$grp" == "GroupID" ]] && continue

  gdir="${groups_dir}/${grp}"
  ref_fa="${gdir}/findref/ref.fa"
  [[ -s "$ref_fa" ]] || { echo "ERROR: missing group ref.fa for ${grp}: $ref_fa" >&2; exit 1; }

  snippy_root="${gdir}/snippy"
  core_dir="${gdir}/core"
  summary_dir="${gdir}/summary"
  mkdir -p "$core_dir" "$summary_dir"

  done_flag="${ROUND2_OUTDIR}/manifest/${grp}.core.done"
  if [[ -f "$done_flag" ]]; then
    echo "[r2_step8] ${grp} core already done -> $done_flag"
    continue
  fi

  # snippy dirs list
  mapfile -t dirs < <(find "$snippy_root" -maxdepth 1 -mindepth 1 -type d | sort)
  [[ "${#dirs[@]}" -gt 0 ]] || { echo "ERROR: no snippy isolate dirs found for ${grp} in $snippy_root" >&2; exit 1; }

  echo "[r2_step8] group=${grp} n=${nmem} dirs=${#dirs[@]}"

  # run snippy-core in core_dir
  cd "$core_dir"
  snippy-core --ref "$ref_fa" "${dirs[@]}" &> "${core_dir}/snippy_core.log"

  # phylo alignment
  snp-sites -b -c -o "${core_dir}/phylo.aln" "${core_dir}/core.full.aln"

  # snp distances -> group summary
  snp-dists "${core_dir}/core.aln" > "${summary_dir}/core.snp-dists.txt"

  # final clustering within group at 50 SNP
  python3 "$CLUSTER_SCRIPT" \
    --snp-dists "${summary_dir}/core.snp-dists.txt" \
    --outdir "$ROUND1_OUTDIR" \
    --threshold 50 \
    --output "${summary_dir}/final_groups.tsv"

  touch "$done_flag"
done < "$groups_tsv"

# -----------------------------
# Merge into FINAL_OUTDIR
# -----------------------------
echo "[r2_step8] merging final outputs into: $FINAL_OUTDIR"

python3 - <<'PY' \
  "${ROUND1_OUTDIR}/summary/initial_groups.tsv" \
  "$groups_tsv" "$groups_dir" "$singletons_tsv" "$FINAL_OUTDIR"
import csv, sys, os
from collections import defaultdict

initial_groups_tsv, groups_tsv, groups_dir, singletons_tsv, final_outdir = sys.argv[1:6]

# isolate -> round1_group
r1_group = {}
with open(initial_groups_tsv, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        r1_group[row["Isolate"].strip()] = row["Group"].strip()

# collect all round2 rows
final_rows = []  # (isolate, species, round2_group)
membership_rows = []  # isolate, round1_group, round2_group

# read each group final_groups.tsv and namespace the within-group labels
with open(groups_tsv) as f:
    next(f)  # header
    for line in f:
        grp, _n = line.strip().split("\t")
        fg = os.path.join(groups_dir, grp, "summary", "final_groups.tsv")
        if not os.path.exists(fg):
            raise SystemExit(f"ERROR: missing {fg}")
        with open(fg, newline="") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for row in r:
                iso = row["Isolate"].strip()
                sp  = row["Species"].strip()
                within = row["Group"].strip()
                r2 = f"{grp}.{within}"
                final_rows.append((iso, sp, r2))
                membership_rows.append((iso, r1_group.get(iso, "NA"), r2))

# add singletons (no round2 processing)
# singletons.tsv: Isolate  Round1Group
singletons = []
if os.path.exists(singletons_tsv) and os.path.getsize(singletons_tsv) > 0:
    with open(singletons_tsv, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            iso = row["Isolate"].strip()
            rg  = row["Round1Group"].strip()
            singletons.append((iso, rg))

# species lookup for singletons: use mash screen top hit from ROUND1_OUTDIR (same logic as cluster script)
def species_from_top3(round1_outdir: str, isolate: str) -> str:
    path = os.path.join(round1_outdir, "mash", isolate, f"{isolate}.mash-screen.top3_hits.txt")
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return "NA"
    with open(path) as f:
        line = f.readline().strip()
    if not line:
        return "NA"
    toks = line.split()
    if len(toks) < 6:
        return "NA"
    desc = " ".join(toks[5:])
    parts = desc.split()
    if len(parts) < 2:
        return "NA"
    return f"{parts[0]} {parts[1]}"

round1_outdir = os.path.dirname(os.path.dirname(initial_groups_tsv))  # .../round1
for iso, rg in singletons:
    sp = species_from_top3(round1_outdir, iso)
    final_rows.append((iso, sp, "Unrelated"))
    membership_rows.append((iso, rg, "Unrelated"))

# sort outputs: Unrelated at bottom
def group_sort_key(g: str):
    if g == "Unrelated":
        return (1, 0, g)
    # sort by initial group number then within-group label
    # g looks like GroupX.GroupY
    try:
        left, right = g.split(".", 1)
    except ValueError:
        left, right = g, ""
    def gnum(x):
        if x.startswith("Group"):
            try: return int(x.replace("Group",""))
            except: return 10**9
        return 10**9
    return (0, gnum(left), gnum(right) if right.startswith("Group") else 10**9, g)

final_rows.sort(key=lambda x: (group_sort_key(x[2]), x[0]))
membership_rows.sort(key=lambda x: (x[1], x[0]))

os.makedirs(final_outdir, exist_ok=True)

with open(os.path.join(final_outdir, "final_groups_all.tsv"), "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["Isolate", "Species", "Group"])
    for iso, sp, g in final_rows:
        w.writerow([iso, sp, g])

with open(os.path.join(final_outdir, "group_membership_map.tsv"), "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["Isolate", "Round1Group", "Round2Group"])
    for iso, g1, g2 in membership_rows:
        w.writerow([iso, g1, g2])
PY

echo "[r2_step8] wrote:"
echo "  ${FINAL_OUTDIR}/final_groups_all.tsv"
echo "  ${FINAL_OUTDIR}/group_membership_map.tsv"
echo "[r2_step8] done"
