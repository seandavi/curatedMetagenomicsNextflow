# Changelog

All notable changes to this pipeline are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
The version is the git tag, the `manifest.version` in `nextflow.config`, and the
workflow revision the orchestrator dispatches — keep all three in lockstep.

## [2.2.1] - 2026-07-04

### Changed
- **Read acquisition is now ENA-first with an SRA fallback** (`fasterq_dump`).
  A class of SRA runs is archived without a QUALITY column, which makes
  `fasterq-dump` exit 3 (`the input data is missing the QUALITY-column`);
  under the retry/`ignore` policy those samples were dropped and dead-lettered.
  All 14 DLQ jobs in the 2.0.7 batch failed this way, yet the reads are
  published and downloadable as FASTQ from EBI/ENA. The process now resolves
  `fastq_ftp` via the ENA `filereport` API, downloads over HTTPS and verifies
  `fastq_md5`, and falls back per run to the original `curl .sra + fasterq-dump`
  path when ENA serves no FASTQ (ingestion lag, submitted-only, controlled
  access). No container change; the downstream output contract is unchanged.
  Note: ENA FASTQ for quality-less runs carries synthetic qualities, so the
  quality-trim step is a no-op for those samples (curation caveat). See ADR-0014.
  Output-neutral in practice: only affects which source serves the bytes, not
  the profiling contract.

### Fixed
- **Resistome KMA failed on 100% of runs (`Error: 2 (No such file or
  directory)`).** KMA writes scratch files under `$TMPDIR`, and the SLURM submit
  templates export `TMPDIR` to a per-job directory that is a sibling of — and not
  bind-mounted alongside — the Nextflow `work` dir. Nextflow propagated that
  `TMPDIR` into the container via `SINGULARITYENV_TMPDIR`, so KMA loaded the CARD
  DB and read the input, then died the moment it touched `$TMPDIR`. Kraken2 and
  MetaPhlAn share the same store/binds but don't use `$TMPDIR`, so only KMA
  tripped. `resistome_kma` now sets `TMPDIR="$PWD"` (the always-mounted task work
  dir).

## [2.2.0] - 2026-07-04

### Removed
- **In-pipeline GTDB conversion.** The SGB→GTDB translation is a static
  relational mapping over the already-published MetaPhlAn profiles and needs no
  per-run compute context, so it is now done as post-processing rather than in
  the pipeline. Removed the `gtdb.nf` module (`metaphlan_to_gtdb`), the
  `sgb_to_gtdb_db` database process, the vendored `bin/cmgd_sgb_to_gtdb.py`, the
  `skip_gtdb` / `sgb2gtdb_url` parameters, and the `skip_gtdb` field in
  `manifest.json`. The per-branch `gtdb/` output subdirectory is no longer
  produced. See [ADR-0013](docs/adr/0013-remove-gtdb-conversion.md) (supersedes
  ADR-0004).

## [2.1.0] - 2026-07-04

### Added
- **Resistome profiling with KMA against CARD**, replacing the RGI/CARD step.
  KMA maps the host-decontaminated reads against a KMA-indexed CARD reference
  (the broadstreet homolog-model FASTA, indexed once into `store_dir` by the new
  `card_kma_db` process) and is cheap enough to run on **both the full and
  rarefied branches** (RGI ran full-only). Outputs are now `card_kma.*.gz`
  (`res`, `mapstat`, `aln`, `fsa`, `frag`) under `resistome/`, replacing the
  `rgi_bwt.*_mapping_data.txt.gz` tables. `rgi_aligner` is removed and
  `card_db_url` now points at the broadstreet tarball. Validated end-to-end
  against real CARD data. See [ADR-0012](docs/adr/0012-resistome-kma-card.md)
  (supersedes ADR-0007).
- **Minimal CI** (`.github/workflows/ci.yml`): a config check across all
  profiles under the current parser plus a container-free `-stub-run`, on every
  push and pull request.

### Changed
- **`nextflow.config` now parses under the strict (v2) Nextflow config parser.**
  The run-lifecycle telemetry hooks were reworked off top-level `import`/`def`
  helpers (which the v2 parser rejects) into inlined `workflow.onComplete` /
  `onError` closures. The pipeline scripts still target the v1 script grammar, so
  runs pin `NXF_SYNTAX_PARSER=v1` (set in the `justfile` recipes and the CI stub
  step) until that separate migration is done.

### Fixed
- **Telemetry `onComplete`/`onError` POSTs now actually fire.** They had silently
  never worked: `ProcessBuilder` copies its arg list into a `String[]`, which
  threw an arraycopy type-mismatch on the interpolated (GString) payload/URL and
  was swallowed by the surrounding try/catch. Fixed by coercing the arg list with
  `*.toString()`. (Pre-existing bug, unrelated to the parser rework.)
## [2.0.7] - 2026-07-02

### Changed
- **Failure isolation: a single bad sample no longer poisons its batch.** The
  default `errorStrategy` for per-sample compute processes now ends in `ignore`
  instead of `finish`. `finish` halted submission of NEW tasks on any non-OOM
  failure, so when one sample's early `fasterq_dump` failed (e.g. exit 3 on a
  bad/unavailable SRA accession) the whole batch's downstream tasks never
  launched — none of the batch-mates reached `MARK_COMPLETE` and all up-to-25
  samples were dead-lettered. Live incident: 4 bad downloads took down 56
  samples, 52 of them healthy collateral. New policy:
  `{ (task.exitStatus in 137..140 && task.attempt <= 4) ? 'retry' : (task.attempt <= 2 ? 'retry' : 'ignore') }`
  — retry the OOM/scheduler-kill family up to 4× (with escalating memory), retry
  any other failure once, then `ignore` (drop just that sample).
- Applied the same retry-then-`ignore` shape to the `google` (Google Batch)
  profile, which previously used `terminate` — strictly worse than `finish`. The
  exit-14 spot/preemption retry is preserved.

### Fixed
- **Shared/critical processes stay fail-hard.** Added `withLabel` fail-hard
  overrides so failure isolation never silently ruins a run:
  - `db_setup` — every reference/database process in
    `modules/processes/databases.nf` (install_metaphlan_db, chocophlan_db,
    utility_mapping_db, uniref_db, sgb_to_gtdb_db, kraken_db, card_db,
    kneaddata_human_database, kneaddata_mouse_database). A broken shared DB must
    stop the run loudly, not silently ruin every sample.
  - `finalize` — MARK_COMPLETE (completion sentinel) and sample_manifest
    (publish/manifest). `ignore`ing these would mark a sample complete without
    outputs or drop its manifest.
  These retry transient failures, then `finish` (never `ignore`).

### Notes
- Profile resolution: the production daemon runs `-profile alpine,gcs`. Neither
  `alpine` nor `gcs` overrides `errorStrategy`, so the effective strategy for
  that run is the `conf/base.config` default — which is exactly what this change
  fixes. The `google` profile's `terminate` override only activates under
  `-profile google` and was fixed for completeness.

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

[2.2.1]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.2.0...2.2.1
[2.0.7]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.0.6...2.0.7
[2.0.6]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.0.5...2.0.6
[2.0.5]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.0.4...2.0.5
[2.0.4]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.0.3...2.0.4
[2.0.3]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.0.2...2.0.3
[2.0.2]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.0.1...2.0.2
[2.0.1]: https://github.com/seandavi/curatedMetagenomicsNextflow/compare/2.0.0...2.0.1
[2.0.0]: https://github.com/seandavi/curatedMetagenomicsNextflow/releases/tag/2.0.0
