# Changelog

## v0.2.0 - 2026-02-25

### Added
- Two-round variant calling workflow:
  - Round 1: cohort-wide SNV calling + initial clustering (default 200 SNP).
  - Round 2: per-initial-group reference selection + SNV calling + final clustering (default 50 SNP).
- CLI options to override clustering thresholds:
  - `--round1-threshold`
  - `--round2-threshold`

### Changed
- Output layout organized into QC/Mash plus two variant-calling rounds, with user-friendly `final_output/`.
- Final outputs emphasize `final_groups_all.tsv` and per-group SNP distance matrices.

### Notes
- Pipeline remains modules-based (Rivanna); containerization is planned for a later release along with instructions for downloading and using database.

## v0.1.0
- First fully operational SLURM-based SNV pipeline
- Mash reference selection implemented
- snippy + snippy-core integration
- SNP distance matrix generation
- Modules-only execution model
