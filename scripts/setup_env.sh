#!/usr/bin/env bash
# Source from each SLURM step. CONDA_ENV is exported by the launcher.
[[ -n "${CONDA_ENV:-}" && -d "${CONDA_ENV}/bin" ]] || { echo "ERROR: CONDA_ENV is unset or invalid: ${CONDA_ENV:-<unset>}" >&2; exit 1; }
export PATH="${CONDA_ENV}/bin:${PATH}"
export CONDA_PREFIX="${CONDA_ENV}"
for tool in python snippy mash; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool unavailable in CONDA_ENV" >&2; exit 1; }
done
