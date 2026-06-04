# Changelog

All notable changes to this pipeline are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
The version is the git tag, the `manifest.version` in `nextflow.config`, and the
workflow revision the orchestrator dispatches — keep all three in lockstep.

## [2.0.6] - 2026-06-04

### Fixed
- `metaphlan_to_gtdb` failed with `unrecognized arguments: -m`. The container
  image (`metaphlan4.2.2`) bakes an **older** `/usr/local/bin/sgb_to_gtdb_profile.py`
  (only `-i`/`-o`) that shadowed the pipeline's vendored `bin/` copy on PATH.
  Renamed the vendored script to `cmgd_sgb_to_gtdb.py` (no name collision) and
  updated the call. (Container should also stop baking the script; the rename is
  the immediate fix.)

## [2.0.5] - 2026-06-04

### Fixed
- `metaphlan_to_gtdb` no longer wrongly reports "no SGB2GTDB mapping table
  found". It now references the downloaded table directly by its known name
  (the basename of `sgb2gtdb_url`) instead of `find`-ing the staged db dir —
  which Nextflow stages as a symlink that `find` won't descend, so the lookup
  found nothing. The new `-s` check also catches an empty/failed download.
  (Latent since the GTDB step was added; first hit now that 2.0.4 runs reach
  this stage.)

## [2.0.4] - 2026-06-04

### Changed
- `errorStrategy` for non-retryable failures is now `finish` instead of
  `terminate`, so a single unrecoverable task lets already-running siblings
  complete rather than aborting the whole batch (ADR-0010).
- Split the GCS storage profile: `gcs` now publishes outputs to GCS only
  (HPC-safe, `-profile alpine,gcs`, local scratch workDir); the new `gcswork`
  profile carries the GCS workDir for cloud compute only (ADR-0011).
- `manifest.version` bumped to `2.0.4` so published paths and telemetry match
  the git tag.

### Fixed
- `google.project` corrected to `curatedmetagenomicdata` (was `omicidx-338300`).

## [2.0.3] - 2026-06-04

### Fixed
- RGI container tag corrected to `6.0.5--pyh05cac1d_0` (was
  `6.0.5--pyha8f3691_0`, which did not exist → `manifest unknown` → the
  `resistome` process failed at image pull and terminated the whole run).

### Documentation
- ADR-0009: retry and failure-handling policy.

## [2.0.2] - 2026-06-04

### Changed
- `rarefy_fastq` now sets explicit resources (`cpus = 2`,
  `memory = { 8.GB * task.attempt }`); previously it ran at the cluster's
  default memory.

## [2.0.1] - 2026-06-04

### Added
- Per-task `.command.out` (stdout) upload via the telemetry `afterScript`, so
  processes that report on stdout (e.g. kraken2) are no longer log-less.

### Changed
- Dockerfile synced with the metaphlan4.2.2 image recipe.

## [2.0.0] - 2026-06-03

Baseline of the 2.x line. Core metagenomic pipeline, with decisions recorded in
`docs/adr/`:

### Added
- Single base container image plus per-tool biocontainers for new tools
  (ADR-0001).
- Dual-branch profiling: full depth + rarefied (ADR-0003).
- GTDB taxonomy conversion via a vendored, `store_dir`-backed mapping table
  (ADR-0004).
- Per-sample provenance & read-accounting manifest (ADR-0005).
- Complementary read-based profiling with Kraken2 + Bracken (ADR-0006).
- Resistome profiling with RGI against CARD (ADR-0007).
- Per-sample post-decontamination FastQC (ADR-0008).

### Notes
- HUMAnN functional profiling is deferred pending MetaPhlAn/HUMAnN version
  alignment (ADR-0002).

[2.0.6]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.0.5...2.0.6
[2.0.5]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.0.4...2.0.5
[2.0.4]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.0.3...2.0.4
[2.0.3]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.0.2...2.0.3
[2.0.2]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.0.1...2.0.2
[2.0.1]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.0.0...2.0.1
[2.0.0]: https://github.com/seandavi/curatedMetagenomicsNextflow/releases/tag/2.0.0
