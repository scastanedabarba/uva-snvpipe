#!/usr/bin/env python3
"""Calculate alignment-wide callable-core statistics from Snippy core.full.aln.

Positions containing -, N, X, or n in any non-reference isolate are excluded.
Reference is intentionally excluded from this mask to match the original script.
"""
import argparse
import csv
from pathlib import Path
from Bio import SeqIO


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('alignment', help='Snippy core.full.aln FASTA')
    p.add_argument('--output', required=True, help='Output TSV file')
    args = p.parse_args()

    records = list(SeqIO.parse(args.alignment, 'fasta'))
    if not records:
        p.error('Alignment is empty')
    references = [r for r in records if r.id == 'Reference']
    if len(references) != 1:
        p.error('Alignment must contain exactly one Reference sequence')
    isolates = [r for r in records if r.id != 'Reference']
    if not isolates:
        p.error('Alignment contains no isolate sequences')
    total = len(references[0].seq)
    if any(len(r.seq) != total for r in records):
        p.error('Alignment sequences have different lengths')

    # Memory-efficient: no per-position pandas DataFrame.
    invalid = set('-NXnx')
    callable_count = sum(all(base not in invalid for base in bases)
                         for bases in zip(*(str(r.seq) for r in isolates)))
    excluded = total - callable_count
    pct = 100.0 * callable_count / total if total else 0.0
    with open(args.output, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['Alignment_Path', 'Total_Positions', 'Callable_Positions',
                    'Excluded_Positions', 'Percent_Callable'])
        w.writerow([str(Path(args.alignment).resolve()), total, callable_count,
                    excluded, f'{pct:.2f}'])
    print(f'[corealn_metrics] {args.output}: {callable_count}/{total} ({pct:.2f}%) callable')


if __name__ == '__main__':
    main()

