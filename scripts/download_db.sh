#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd -P)"
CONFIG=""
MANIFEST="$REPO_DIR/databases/db_manifest.tsv"

DB_NAME="uva_eskape_refseq_db"
DB_VERSION=""
OUTDIR=""
FORCE=0
LOCAL=0

usage() {
  cat <<'USAGE'
Usage: bash scripts/download_db.sh --outdir PATH [--config FILE] [--name NAME] [--version VERSION] [--force] [--local]

  --outdir PATH      Required parent directory for databases (installed under PATH/databases/)
  --config FILE      Configuration file (default: repo/config/config.sh)
  --name NAME        Database name in manifest (default: uva_eskape_refseq_db)
  --version VERSION  Database version in manifest (default: first matching entry)
  --force            Replace an existing installation
  --local            Run directly instead of submitting a SLURM job
  -h, --help         Show this help

Example:
  bash scripts/download_db.sh --outdir /scratch/username/snvpipeline
USAGE
}
die() { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"; }

while (($#)); do
  case "$1" in
    --config|--outdir|--name|--version)
      (($# >= 2)) || die "Missing value for $1"
      [[ -n "$2" ]] || die "Empty value for $1"
      case "$1" in
        --config) CONFIG="$2" ;;
        --outdir) OUTDIR="$2" ;;
        --name) DB_NAME="$2" ;;
        --version) DB_VERSION="$2" ;;
      esac
      shift 2 ;;
    --force) FORCE=1; shift ;;
    --local) LOCAL=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

# Capture command-line values before sourcing site configuration.
CLI_OUTDIR="$OUTDIR"
[[ -n "$CLI_OUTDIR" ]] || die "Missing required --outdir PATH"
CONFIG="${CONFIG:-$REPO_DIR/config/config.sh}"
[[ -f "$CONFIG" && -r "$CONFIG" ]] || die "Configuration file not found or unreadable: $CONFIG"
CONFIG="$(cd "$(dirname "$CONFIG")" && pwd -P)/$(basename "$CONFIG")"
# shellcheck source=/dev/null
source "$CONFIG"
OUTDIR="$CLI_OUTDIR"

[[ -f "$MANIFEST" ]] || die "Manifest not found: $MANIFEST"
need awk
need sha256sum
need tar
if [[ -n "${CONDA_ENV:-}" && -x "${CONDA_ENV}/bin/aria2c" ]]; then
  export PATH="${CONDA_ENV}/bin:$PATH"
fi
if command -v aria2c >/dev/null 2>&1; then
  DOWNLOADER=aria2
elif command -v curl >/dev/null 2>&1; then
  DOWNLOADER=curl
elif command -v wget >/dev/null 2>&1; then
  DOWNLOADER=wget
else
  die "Requires aria2c, curl, or wget"
fi

# Normalize parent without requiring it to exist; --outdir is never the final DB directory.
mkdir -p "$(dirname "$OUTDIR")"
OUTDIR="$(cd "$(dirname "$OUTDIR")" && pwd -P)/$(basename "$OUTDIR")"
[[ "$OUTDIR" != / ]] || die "Refusing to install into /"
BASE_DIR="$OUTDIR"

row="$(awk -v name="$DB_NAME" -v ver="$DB_VERSION" 'NR>1 && $1==name && (ver=="" || $2==ver) {print; exit}' "$MANIFEST")"
[[ -n "$row" ]] || die "No manifest entry for name='$DB_NAME' version='${DB_VERSION:-any}'"
read -r found_name found_version url expected_sha _ <<< "$row"
[[ -n "$url" && -n "$expected_sha" ]] || die "Manifest row lacks URL or SHA256"
[[ "$expected_sha" =~ ^[[:xdigit:]]{64}$ ]] || die "Invalid SHA256 in manifest"
archive_name="$(basename "${url%%\?*}")"
[[ "$archive_name" == *.tar.gz ]] || die "Expected .tar.gz URL, got $archive_name"
archive_basename="${archive_name%.tar.gz}"
[[ "$archive_basename" =~ ^[a-zA-Z0-9][a-zA-Z0-9._-]*$ ]] || die "Unsafe archive basename: $archive_basename"
OUTDIR="$BASE_DIR/databases/$archive_basename"
mkdir -p "$BASE_DIR/databases"

if (( ! LOCAL )); then
  need sbatch
  LOGDIR="$BASE_DIR/databases/logs"
  mkdir -p "$LOGDIR"
  args=(--local --config "$CONFIG" --outdir "$BASE_DIR" --name "$DB_NAME")
  [[ -z "$DB_VERSION" ]] || args+=(--version "$DB_VERSION")
  (( ! FORCE )) || args+=(--force)
  printf -v job_cmd '%q ' bash "$SCRIPT_DIR/download_db.sh" "${args[@]}"
  sbatch_args=(--parsable --job-name=snv-db-download --output="$LOGDIR/download-%j.out" --error="$LOGDIR/download-%j.err")
  [[ -z "${ACCOUNT:-}" ]] || sbatch_args+=(--account="$ACCOUNT")
  [[ -z "${PARTITION:-}" ]] || sbatch_args+=(--partition="$PARTITION")
  sbatch_args+=(--time="${DB_DOWNLOAD_TIME:-24:00:00}" --cpus-per-task=2 --mem=8G)
  job_id="$(sbatch "${sbatch_args[@]}" --wrap="$job_cmd")"
  echo "[download_db] Submitted SLURM job: $job_id"
  echo "[download_db] Logs: $LOGDIR/download-${job_id}.out and .err"
  echo "[download_db] Track with: squeue -j $job_id"
  exit 0
fi

echo "[download_db] Installation target: $OUTDIR"
if [[ -d "$OUTDIR" ]] && [[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  if (( ! FORCE )); then
    if [[ -d "$OUTDIR/mash" ]]; then
      echo "[download_db] Database directory already exists: $OUTDIR"
      echo "[download_db] Skipping download; use --force to reinstall."
      exit 0
    fi
    die "Database target is nonempty but missing mash/: $OUTDIR (inspect it before replacing)"
  fi
fi

DOWNLOAD_DIR="$BASE_DIR/databases/.downloads"
mkdir -p "$DOWNLOAD_DIR"
archive="$DOWNLOAD_DIR/$archive_name"
TMPDIR_INSTALL="$(mktemp -d "$(dirname "$OUTDIR")/.uva-snvdb.XXXXXXXX")"
trap 'rm -rf -- "$TMPDIR_INSTALL"' EXIT
staging="$TMPDIR_INSTALL/extracted"
mkdir -p "$staging"

archive_valid=0
if [[ -f "$archive" ]]; then
  echo "[download_db] Checking existing archive: $archive"
  actual_sha="$(sha256sum "$archive" | awk '{print $1}')"
  if [[ "${actual_sha,,}" == "${expected_sha,,}" ]]; then
    echo "[download_db] Existing archive SHA256 OK; skipping download."
    archive_valid=1
  else
    echo "[download_db] Existing archive is incomplete or invalid; resuming download."
  fi
fi

if (( ! archive_valid )); then
  echo "[download_db] Downloading $found_name version $found_version"
  if [[ "$DOWNLOADER" == aria2 ]]; then
    aria2c --continue=true --max-connection-per-server=8 --split=8 \
      --min-split-size=4M --max-tries=30 --retry-wait=15 \
      --timeout=120 --connect-timeout=30 --file-allocation=none \
      --dir="$DOWNLOAD_DIR" --out="$archive_name" "$url"
  elif [[ "$DOWNLOADER" == curl ]]; then
    curl --fail --location --continue-at - --retry 30 --retry-delay 15 \
      --retry-all-errors --output "$archive" "$url"
  else
    wget --continue --tries=30 --waitretry=15 --output-document="$archive" "$url"
  fi
  echo "[download_db] Verifying SHA256..."
  actual_sha="$(sha256sum "$archive" | awk '{print $1}')"
  [[ "${actual_sha,,}" == "${expected_sha,,}" ]] || die "SHA256 mismatch: expected $expected_sha, got $actual_sha"
  echo "[download_db] SHA256 OK"
fi

# Validate archive before extraction. Reject links, special files, absolute paths,
# traversal, and members outside the expected top-level directory.
python_cmd="$(command -v python3 || command -v python || true)"
[[ -n "$python_cmd" ]] || die "Python is required to validate the archive"
"$python_cmd" - "$archive" "$archive_basename" <<'PYVALIDATE'
import sys
import tarfile
from pathlib import PurePosixPath

archive, expected_root = sys.argv[1:]
seen = False
with tarfile.open(archive, 'r:gz') as tf:
    for member in tf:
        name = member.name
        normalized = name.rstrip('/')
        parts = PurePosixPath(normalized).parts
        if (not normalized or name.startswith('/') or
                '..' in parts or '.' in parts or
                parts[0] != expected_root or
                not (member.isfile() or member.isdir()) or
                (len(parts) == 1 and not member.isdir())):
            raise SystemExit(f'Unsafe or unexpected archive member: {name}')
        seen = True
if not seen:
    raise SystemExit('Archive contains no members')
PYVALIDATE

echo "[download_db] Extracting database..."
tar -xzf "$archive" -C "$staging" --strip-components=1
[[ -d "$staging/mash" ]] || die "Extracted database missing mash/ directory"

if [[ -e "$OUTDIR" ]]; then
  if (( FORCE )); then
    echo "[download_db] Replacing existing installation: $OUTDIR"
    rm -rf -- "$OUTDIR"
  else
    rmdir "$OUTDIR" 2>/dev/null || die "Destination unexpectedly became nonempty: $OUTDIR"
  fi
fi
mv "$staging" "$OUTDIR"
rm -f -- "$archive" "$archive.aria2"
echo "[download_db] Installed: $OUTDIR"
echo "[download_db] Name: $found_name | Version: $found_version"
echo "[download_db] Set DB_DIR=\"$OUTDIR\" in config/config.sh to use this installation."

