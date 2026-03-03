#!/usr/bin/env bash
set -euo pipefail

MANIFEST="databases/db_manifest.tsv"

# Defaults
DB_NAME="uva_eskape_refseq_db"
DB_VERSION=""                  # if empty, first match by name
DB_ROOT="databases"            # install root (default: repo-local ./databases)
FORCE=0                        # set to 1 to overwrite
DELETE_ARCHIVE=1   # default: delete tar after successful install

usage() {
  cat <<'EOF'
Usage:
  bash scripts/download_db.sh [--name DB_NAME] [--version DB_VERSION] [--db-root PATH] [--force]

Defaults:
  --name      uva_eskape_refseq_db
  --version   (first matching version in manifest)
  --db-root   databases
  --force     false

Examples:
  # Default: installs into ./databases/<archive_basename> (e.g., ./databases/uva_eskape_2026-01-23)
  bash scripts/download_db.sh

  # Install a specific version
  bash scripts/download_db.sh --version 2026-01-23

  # Install somewhere else (e.g., scratch)
  bash scripts/download_db.sh --db-root /scratch/$USER/snvdb

  # Reinstall (overwrite)
  bash scripts/download_db.sh --force
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"; }

need awk
need sha256sum
need tar

DOWNLOADER=""
if command -v curl >/dev/null 2>&1; then
  DOWNLOADER="curl"
elif command -v wget >/dev/null 2>&1; then
  DOWNLOADER="wget"
else
  die "Need curl or wget"
fi

# -------------------------
# Parse args
# -------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --name) DB_NAME="${2:-}"; shift 2;;
    --version) DB_VERSION="${2:-}"; shift 2;;
    --db-root) DB_ROOT="${2:-}"; shift 2;;
    --force) FORCE=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) die "Unknown arg: $1";;
  esac
done

[[ -f "$MANIFEST" ]] || die "Manifest not found: $MANIFEST"

# -------------------------
# Read manifest row
# -------------------------
row="$(
  awk -v name="$DB_NAME" -v ver="$DB_VERSION" '
    NR==1 {next}
    $1==name && (ver=="" || $2==ver) {print $0; exit 0}
  ' "$MANIFEST"
)"

db_name="$(echo "$row" | awk '{print $1}')"
db_version="$(echo "$row" | awk '{print $2}')"
url="$(echo "$row" | awk '{print $3}')"
sha256_expected="$(echo "$row" | awk '{print $4}')"

[[ -n "$url" ]] || die "Empty URL in manifest"
[[ -n "$sha256_expected" ]] || die "Empty sha256 in manifest"

# Determine archive filename from URL (Zenodo uses .../files/<FILENAME>?download=1)
url_file="$(basename "${url%%\?*}")"                 # strip querystring, keep filename
archive_basename="${url_file%.tar.gz}"               # e.g., uva_eskape_2026-01-23
[[ "$archive_basename" != "$url_file" ]] || die "Expected a .tar.gz URL, got: $url_file"

# Install target (clean, non-nested)
dest_dir="${DB_ROOT}/${archive_basename}"
archive_path="${DB_ROOT}/${url_file}"

mkdir -p "$DB_ROOT"

# -------------------------
# Download
# -------------------------
if [[ -f "$archive_path" && "$FORCE" != "1" ]]; then
  echo "[download_db] Using existing archive: $archive_path (use --force to re-download)"
else
  echo "[download_db] Downloading: $url"
  if [[ "$DOWNLOADER" == "curl" ]]; then
    curl -L --fail -o "$archive_path" "$url"
  else
    wget -O "$archive_path" "$url"
  fi
fi

# -------------------------
# Verify checksum
# -------------------------
echo "[download_db] Verifying sha256..."
sha256_actual="$(sha256sum "$archive_path" | awk '{print $1}')"
[[ "$sha256_actual" == "$sha256_expected" ]] || die "SHA256 mismatch
Expected: $sha256_expected
Actual:   $sha256_actual"

echo "[download_db] sha256 OK"

# -------------------------
# Extract (de-nest)
# -------------------------
if [[ -d "$dest_dir" && "$FORCE" != "1" ]] && find "$dest_dir" -mindepth 1 -maxdepth 1 | read -r _; then
  echo "[download_db] Destination not empty: $dest_dir (use --force to overwrite)"
else
  if [[ "$FORCE" == "1" ]]; then
    echo "[download_db] --force: removing existing $dest_dir"
    rm -rf "$dest_dir"
  fi

  mkdir -p "$dest_dir"

  echo "[download_db] Extracting -> $dest_dir"
  tar -xzf "$archive_path" -C "$dest_dir" --strip-components=1

  # Sanity check: ensure key directory exists
  if [[ ! -d "$dest_dir/mash" ]]; then
    echo "ERROR: Extraction failed — expected directory $dest_dir/mash not found" >&2
    exit 1
  fi

  echo "[download_db] Extraction complete"
fi

# -------------------------
# Remove archive
# -------------------------
echo "[download_db] Removing archive: $archive_path"
rm -f "$archive_path"

echo "[download_db] Database installed at: $dest_dir"
echo "  name:    $db_name"
echo "  version: $db_version"
echo "  path:    $dest_dir"
