#!/usr/bin/env bash
#SBATCH --job-name=r2_core
# Environment is exported explicitly by the launcher

set -euo pipefail

: "${PIPELINE_ROOT:?ERROR: PIPELINE_ROOT is not set}"
SCRIPT_DIR="${PIPELINE_ROOT}/workflow"
source "${SCRIPT_DIR}/../scripts/setup_env.sh"

OUTDIR_ROOT="$(readlink -f "$1")"
ROUND2_OUTDIR="$(readlink -f "$2")"
FINAL_OUTDIR="$(readlink -f "$3")"
CLUSTER_SCRIPT="$(readlink -f "$4")"
THRESHOLD="${5:-50}"
METRICS_SCRIPT="${PIPELINE_ROOT}/scripts/corealn_metrics.py"

groups_tsv="${ROUND2_OUTDIR}/manifest/groups.tsv"
singletons_tsv="${ROUND2_OUTDIR}/manifest/singletons.tsv"
groups_dir="${ROUND2_OUTDIR}/groups"
round2_summary="${ROUND2_OUTDIR}/summary"
round1_summary="${OUTDIR_ROOT}/variant_calling_round1/summary"
final_round1="${FINAL_OUTDIR}/round1"
final_round2="${FINAL_OUTDIR}/round2"

[[ -s "$groups_tsv" ]] || { echo "ERROR: missing groups.tsv: $groups_tsv" >&2; exit 1; }
[[ -f "$CLUSTER_SCRIPT" ]] || { echo "ERROR: missing cluster script: $CLUSTER_SCRIPT" >&2; exit 1; }
[[ -f "$METRICS_SCRIPT" ]] || { echo "ERROR: missing metrics script: $METRICS_SCRIPT" >&2; exit 1; }

mkdir -p "$FINAL_OUTDIR" "$round2_summary" "$final_round1" "$final_round2"

# Preserve round-one matrix and callability report in the final deliverables.
for report in core.snp-dists.txt core_alignment_metrics.tsv; do
  [[ -s "${round1_summary}/${report}" ]] || { echo "ERROR: missing round1 report: ${round1_summary}/${report}" >&2; exit 1; }
  cp -f "${round1_summary}/${report}" "${final_round1}/${report}"
done

echo "[r2_step8] per-group core + final clustering ..."
while read -r grp nmem; do
  [[ "$grp" == "GroupID" || -z "$grp" ]] && continue

  gdir="${groups_dir}/${grp}"
  ref_fa="${gdir}/findref/ref.fa"
  [[ -s "$ref_fa" ]] || { echo "ERROR: missing group ref.fa for ${grp}: $ref_fa" >&2; exit 1; }

  snippy_root="${gdir}/snippy"
  done_flag="${ROUND2_OUTDIR}/manifest/${grp}.core.done"
  group_snp="${round2_summary}/${grp}.core.snp-dists.txt"
  group_metrics="${round2_summary}/${grp}.core_alignment_metrics.tsv"
  group_final="${round2_summary}/${grp}.final_groups.tsv"

  if [[ -f "$done_flag" ]]; then
    echo "[r2_step8] ${grp} core already done -> $done_flag"
    # Backfill metrics for groups completed by an older version of the workflow.
    if [[ ! -s "$group_metrics" ]]; then
      [[ -s "${snippy_root}/core.full.aln" ]] || { echo "ERROR: missing core.full.aln for ${grp}" >&2; exit 1; }
      python3 "$METRICS_SCRIPT" "${snippy_root}/core.full.aln" --output "$group_metrics"
    fi
    [[ -s "$group_snp" && -s "$group_final" ]] || { echo "ERROR: missing completed group reports for ${grp}" >&2; exit 1; }
    cp -f "$group_snp" "$final_round2/$(basename "$group_snp")"
    cp -f "$group_metrics" "$final_round2/$(basename "$group_metrics")"
    continue
  fi

  mapfile -t dirs < <(find "$snippy_root" -maxdepth 1 -mindepth 1 -type d | sort)
  [[ "${#dirs[@]}" -gt 0 ]] || { echo "ERROR: no snippy isolate dirs found for ${grp} in $snippy_root" >&2; exit 1; }

  echo "[r2_step8] group=${grp} n=${nmem} dirs=${#dirs[@]}"

  cd "$snippy_root"
  snippy-core --ref "$ref_fa" "${dirs[@]}" &> "${snippy_root}/snippy_core.log"

  snp-sites -b -c -o "${snippy_root}/phylo.aln" "${snippy_root}/core.full.aln"

  python3 "$METRICS_SCRIPT" "${snippy_root}/core.full.aln" --output "$group_metrics"

  snp-dists "${snippy_root}/core.aln" > "$group_snp"

  python3 "$CLUSTER_SCRIPT" \
    --snp-dists "$group_snp" \
    --outdir "$OUTDIR_ROOT" \
    --threshold "$THRESHOLD" \
    --output "$group_final"

  cp -f "$group_snp" "$final_round2/$(basename "$group_snp")"
  cp -f "$group_metrics" "$final_round2/$(basename "$group_metrics")"

  touch "$done_flag"
done < "$groups_tsv"

echo "[r2_step8] merging final outputs into: $FINAL_OUTDIR"

python3 - <<'PY' \
  "${OUTDIR_ROOT}/variant_calling_round1/summary/initial_groups.tsv" \
  "$groups_tsv" "$round2_summary" "$singletons_tsv" "$FINAL_OUTDIR" "$OUTDIR_ROOT"
import csv, sys, os

initial_groups_tsv, groups_tsv, round2_summary, singletons_tsv, final_outdir, top_outdir = sys.argv[1:7]

def species_from_top3(outdir: str, isolate: str) -> str:
    path = os.path.join(outdir, "mash", isolate, f"{isolate}.mash-screen.top3_hits.txt")
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

r1_group = {}
with open(initial_groups_tsv, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        r1_group[row["Isolate"].strip()] = row["Group"].strip()

final_rows = []
with open(groups_tsv) as f:
    next(f)
    for line in f:
        grp, _n = line.strip().split("\t")
        fg = os.path.join(round2_summary, f"{grp}.final_groups.tsv")
        if not os.path.exists(fg):
            raise SystemExit(f"ERROR: missing {fg}")
        with open(fg, newline="") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for row in r:
                iso = row["Isolate"].strip()
                sp = row["Species"].strip()
                secondary = row["Group"].strip()
                final_rows.append((iso, sp, grp, secondary))

if os.path.exists(singletons_tsv) and os.path.getsize(singletons_tsv) > 0:
    with open(singletons_tsv, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            iso = row["Isolate"].strip()
            sp = species_from_top3(top_outdir, iso)
            final_rows.append((iso, sp, "Unrelated", "Unrelated"))

def grp_num(g):
    if g.startswith("Group"):
        try: return int(g.replace("Group",""))
        except: return 10**9
    if g == "Unrelated":
        return 10**9
    return 10**9

final_rows.sort(key=lambda x: (grp_num(x[2]), x[2], grp_num(x[3]), x[3], x[0]))

os.makedirs(final_outdir, exist_ok=True)
out_path = os.path.join(final_outdir, "final_groups_all.tsv")
with open(out_path, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["Isolate", "Species", "Primary_Group", "Secondary_Group"])
    for iso, sp, pg, sg in final_rows:
        w.writerow([iso, sp, pg, sg])
PY

echo "[r2_step8] wrote:"
echo "  ${FINAL_OUTDIR}/final_groups_all.tsv"
echo "  ${final_round1}/core.snp-dists.txt"
echo "  ${final_round1}/core_alignment_metrics.tsv"
echo "  ${final_round2}/*.core.snp-dists.txt"
echo "  ${final_round2}/*.core_alignment_metrics.tsv"
echo "[r2_step8] done"

