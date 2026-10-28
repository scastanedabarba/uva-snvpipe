#!/usr/bin/env python3
"""Write analysis software versions and SNP clustering thresholds as a TSV.

Reads installed package versions from the specified Conda prefix's conda-meta,
so no activation is required. Standard-library only.
"""
import argparse
import csv
import json
import subprocess
from pathlib import Path

PACKAGES = (
    'python', 'snippy', 'mash', 'snp-dists', 'snp-sites',
    'trim-galore', 'trimmomatic', 'cutadapt', 'fastqc', 'openjdk',
    'numpy', 'python-dateutil', 'aria2', 'curl',
)


def package_versions(prefix):
    metadata = Path(prefix) / 'conda-meta'
    if not metadata.is_dir():
        raise SystemExit(f'ERROR: Conda metadata directory not found: {metadata}')
    versions = {}
    for file in metadata.glob('*.json'):
        with file.open(encoding='utf-8') as fh:
            record = json.load(fh)
        name, version = record.get('name'), record.get('version')
        if name in PACKAGES and version:
            versions[name] = version
    return versions


def git_value(repo, *args):
    try:
        return subprocess.check_output(
            ['git', '-C', str(repo), *args], stderr=subprocess.DEVNULL,
            text=True, timeout=10,
        ).strip() or 'NA'
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return 'NA'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, help='Output TSV file')
    p.add_argument('--conda-env', required=True, help='Installed Conda prefix')
    p.add_argument('--repo', required=True, help='Pipeline repository root')
    p.add_argument('--db-path', required=True, help='Selected database directory')
    p.add_argument('--round1-threshold', required=True, type=int)
    p.add_argument('--round2-threshold', required=True, type=int)
    args = p.parse_args()
    if args.round1_threshold < 0 or args.round2_threshold < 0:
        p.error('SNP thresholds must be nonnegative')

    versions = package_versions(args.conda_env)
    repo = Path(args.repo).resolve()
    db = Path(args.db_path).resolve()
    rows = [
        ('pipeline', 'git_commit', git_value(repo, 'rev-parse', 'HEAD')),
        ('pipeline', 'git_version', git_value(repo, 'describe', '--tags', '--always')),
        ('database', 'path', str(db)),
        ('database', 'release_directory', db.name),
        ('parameters', 'round1_snp_threshold', str(args.round1_threshold)),
        ('parameters', 'round2_snp_threshold', str(args.round2_threshold)),
        ('environment', 'conda_prefix', str(Path(args.conda_env).resolve())),
    ]
    rows.extend(('software', name, versions.get(name, 'NOT_FOUND')) for name in PACKAGES)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(('Category', 'Name', 'Value'))
        writer.writerows(rows)
    print(f'[provenance] wrote {output}')


if __name__ == '__main__':
    main()

