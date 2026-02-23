#!/usr/bin/env python3
"""
cluster_isolates.py

Cluster isolates into groups where members share <= threshold SNP distance
(based on snp-dists matrix), ignore specified labels (e.g., Reference),
and attach speciation results from Mash outputs.

Preferred species source:
  <outdir>/mash/<isolate>/<isolate>.mash-screen.top3_hits.txt
    (uses top hit description; outputs first two words as "Genus species")

Fallback:
  <outdir>/mash/<isolate>/<isolate>.mash-ref.txt
    (uses description field; outputs first two words as "Genus species")

Output TSV columns:
  Isolate, Species, Group
"""

from __future__ import annotations

import argparse
import os
from collections import deque
from typing import Dict, List, Set, Tuple


def read_snp_dists_matrix(path: str, ignore_names: Set[str]) -> Tuple[List[str], Dict[str, Dict[str, int]]]:
    """
    Reads a snp-dists text matrix. Handles header like:
      snp-dists 0.8.2 A B C Reference
    Returns:
      isolates (filtered, excludes ignore_names)
      dist[row][col] = int distance for filtered names only
    """
    with open(path, "r") as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]

    if not lines:
        raise ValueError(f"Empty SNP dists file: {path}")

    header_tokens = lines[0].split()

    # Handle the common snp-dists header prefix "snp-dists <version>"
    if len(header_tokens) >= 3 and header_tokens[0] == "snp-dists":
        colnames = header_tokens[2:]
    else:
        colnames = header_tokens

    raw_rows: Dict[str, Dict[str, int]] = {}
    for ln in lines[1:]:
        toks = ln.split()
        row = toks[0]
        vals = toks[1:]
        if len(vals) != len(colnames):
            raise ValueError(
                f"Row length mismatch for {row}: got {len(vals)} values, expected {len(colnames)}.\n"
                f"Line: {ln}"
            )
        raw_rows[row] = {colnames[i]: int(vals[i]) for i in range(len(colnames))}

    # Names present as both rows and columns
    names = [n for n in colnames if n in raw_rows]
    filtered = [n for n in names if n not in ignore_names]

    dist: Dict[str, Dict[str, int]] = {}
    for r in filtered:
        dist[r] = {}
        for c in filtered:
            dist[r][c] = raw_rows[r][c]

    return filtered, dist


def build_graph(isolates: List[str], dist: Dict[str, Dict[str, int]], threshold: int) -> Dict[str, Set[str]]:
    """
    Create adjacency list graph where edge exists if dist[i][j] <= threshold and i != j
    """
    adj: Dict[str, Set[str]] = {i: set() for i in isolates}
    n = len(isolates)
    for i in range(n):
        a = isolates[i]
        for j in range(i + 1, n):
            b = isolates[j]
            if dist[a][b] <= threshold:
                adj[a].add(b)
                adj[b].add(a)
    return adj


def connected_components(adj: Dict[str, Set[str]]) -> List[Set[str]]:
    """
    Return connected components of an undirected graph adjacency list.
    """
    seen: Set[str] = set()
    comps: List[Set[str]] = []

    for node in adj:
        if node in seen:
            continue
        q = deque([node])
        seen.add(node)
        comp = {node}
        while q:
            cur = q.popleft()
            for nb in adj[cur]:
                if nb not in seen:
                    seen.add(nb)
                    q.append(nb)
                    comp.add(nb)
        comps.append(comp)

    # Sort: largest first, then lexicographically
    comps.sort(key=lambda s: (-len(s), sorted(s)[0]))
    return comps


def genus_species_from_description(desc: str) -> str:
    """
    Convert a Mash description string into "Genus species" by taking first two words.
    Example:
      "Klebsiella pneumoniae strain CAV1193 chromosome, complete genome"
        -> "Klebsiella pneumoniae"
    """
    desc = desc.strip()
    if not desc:
        return "NA"
    parts = desc.split()
    if len(parts) < 2:
        return "NA"
    return f"{parts[0]} {parts[1]}"


def load_species_from_top3(outdir: str, isolate: str) -> str:
    """
    Reads species from the FIRST LINE of:
      <outdir>/mash/<isolate>/<isolate>.mash-screen.top3_hits.txt

    Expected columns (whitespace-delimited):
      id  shared-hashes  median-multiplicity  p-value  accession|contig  description...

    Returns: "Genus species" or "NA"
    """
    top3 = os.path.join(outdir, "mash", isolate, f"{isolate}.mash-screen.top3_hits.txt")
    if not os.path.isfile(top3) or os.path.getsize(top3) == 0:
        return "NA"

    with open(top3, "r") as fh:
        line = fh.readline().strip()
    if not line:
        return "NA"

    toks = line.split()
    # description starts at token 5 onward (0-based index 5) given your example:
    # 0:identity 1:hashes 2:median_mult 3:pval 4:accession 5+:description
    if len(toks) < 6:
        return "NA"

    desc = " ".join(toks[5:])
    return genus_species_from_description(desc)


def load_species_from_mash_ref(outdir: str, isolate: str) -> str:
    """
    Fallback: reads species from:
      <outdir>/mash/<isolate>/<isolate>.mash-ref.txt

    Your pipeline writes:
      sample <tab> refaccid <tab> refstr...
    We'll take refstr and return first two words.
    """
    mash_ref = os.path.join(outdir, "mash", isolate, f"{isolate}.mash-ref.txt")
    if not os.path.isfile(mash_ref) or os.path.getsize(mash_ref) == 0:
        return "NA"

    with open(mash_ref, "r") as fh:
        line = fh.readline().rstrip("\n")
    if not line:
        return "NA"

    parts = line.split("\t")
    if len(parts) < 3:
        parts = line.split()
    if len(parts) < 3:
        return "NA"

    refstr = " ".join(parts[2:]).strip()
    return genus_species_from_description(refstr)


def load_species(outdir: str, isolate: str) -> str:
    """
    Prefer top3 hits (screen output). Fall back to mash-ref.
    """
    sp = load_species_from_top3(outdir, isolate)
    if sp != "NA":
        return sp
    return load_species_from_mash_ref(outdir, isolate)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Cluster isolates by SNP distance threshold and attach Genus-species calls from Mash outputs."
    )
    ap.add_argument("--snp-dists", required=True, help="Path to core.snp-dists.txt (snp-dists matrix)")
    ap.add_argument("--outdir", required=True, help="Pipeline OUTDIR (contains mash/<isolate>/...)")
    ap.add_argument("--threshold", type=int, default=50, help="SNP threshold for linking isolates (default: 50)")
    ap.add_argument("--ignore", default="Reference", help="Comma-separated names to ignore (default: Reference)")
    ap.add_argument("--output", required=True, help="Output TSV path")
    args = ap.parse_args()

    ignore_names = {x.strip() for x in args.ignore.split(",") if x.strip()}

    isolates, dist = read_snp_dists_matrix(args.snp_dists, ignore_names=ignore_names)
    if not isolates:
        raise SystemExit("ERROR: No isolates found after filtering ignore names.")

    adj = build_graph(isolates, dist, threshold=args.threshold)
    comps = connected_components(adj)

    # Group assignment:
    # - components with size >= 2: Group1, Group2, ...
    # - singletons: Unrelated
    isolate_to_group: Dict[str, str] = {}
    group_idx = 1
    for comp in comps:
        if len(comp) >= 2:
            gid = f"Group{group_idx}"
            group_idx += 1
            for iso in comp:
                isolate_to_group[iso] = gid
        else:
            iso = next(iter(comp))
            isolate_to_group[iso] = "Unrelated"


    # Prepare sorted output rows
    rows = []
    for iso in isolates:
        species = load_species(args.outdir, iso)
        group = isolate_to_group.get(iso, "Unrelated")
        rows.append((iso, species, group))

    def group_sort_key(group_name: str):
        if group_name.startswith("Group"):
            try:
                return (0, int(group_name.replace("Group", "")))
            except ValueError:
                return (0, 999999)
        elif group_name == "Unrelated":
            return (1, 0)
        else:
            return (2, group_name)

    rows.sort(key=lambda x: (group_sort_key(x[2]), x[0]))

    with open(args.output, "w") as out:
        out.write("Isolate\tSpecies\tGroup\n")
        for iso, species, group in rows:
            out.write(f"{iso}\t{species}\t{group}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
