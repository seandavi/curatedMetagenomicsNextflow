# Changelog

All notable changes to this pipeline are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
The version is the git tag, the `manifest.version` in `nextflow.config`, and the
workflow revision the orchestrator dispatches — keep all three in lockstep.

## [Unreleased]

### Fixed
- The HUMAnN-side MetaPhlAn steps (`metaphlan_db_humann`,
  `metaphlan_for_humann`) passed `--db_dir`, which MetaPhlAn 4.1.x rejects
  (`unrecognized arguments: --db_dir`; 4.1 calls it `--bowtie2db`). The option
  is now the bundle field `metaphlan_db_option`. Found by the Alpine pilot
  (#88).
- `humann` passed `--utility-database`, which HUMAnN 3.9 doesn't have
  (`unrecognized arguments`; only 4.0 alpha added it). It is now the bundle
  field `humann_utility_db_option`, null for `humann3.9`. Found by the Alpine
  pilot (#88).
- The recorded HUMAnN version was `SyntaxWarning:`: the 3.9 container prints
  Python 3.12 `SyntaxWarning`s before `humann v3.9`, and the version capture
  took the second word of the merged output. It now reads the `humann` line.

### Added
- **`humann4.0.0a1` bundle** (#87): HUMAnN 4.0.0a1 with MetaPhlAn 4.1.1 and the
  `mpa_vOct22_CHOCOPhlAnSGB_202403` index, plus the bundle's v4-alpha
  ChocoPhlAn/UniRef/utility-mapping databases. HUMAnN 4.0.0a1 has no usable
  upstream image, so `docker/humann4a/` is built and published to
  `ghcr.io/seandavi/humann:4.0.0a1` by `.github/workflows/humann4a-image.yml`
  (ADR-0017). `humann3.9` stays the default; `skip_humann` stays `true`.
- **HUMAnN subworkflow driven by version-pinned bundles** (ADR-0016, #85).
  `--humann_bundle` (default `humann3.9`, registry in
  `conf/humann_bundles.config`) selects the HUMAnN and MetaPhlAn containers,
  the MetaPhlAn version/index and the ChocoPhlAn/UniRef/utility-mapping DB
  names as one unit; an unknown bundle fails at start-up listing the valid
  names (only when `--skip_humann false`). The new `HUMANN` subworkflow runs
  `metaphlan_for_humann` (MetaPhlAn 4.1.1, `mpa_vJun23_CHOCOPhlAnSGB_202307`)
  then `humann` on the full-depth branch only, with its own databases staged
  by `DATABASES` (so `--databases_only --skip_humann false` fetches them).
  Output is published to `<sample>/humann/<bundle>/` with the profile in
  `metaphlan/`; the published MetaPhlAn 4.2.2 taxonomy is unchanged.
  `skip_humann` stays `true` by default. The manifest records the HUMAnN-pass
  versions under `metaphlan_humann`, `bowtie2_humann` and `humann`, plus
  `parameters.humann_bundle` and `humann_metaphlan_index`.
- `humann_maxforks` (default 4) throttles concurrent `humann` tasks.
- **`--databases_only`** pre-stages the reference databases into `store_dir`
  without sample inputs or per-sample processes, so downloads no longer have
  to happen inside a production batch. Database processes are now invoked
  only from a `DATABASES` subworkflow (`modules/subworkflows/databases.nf`).
  See #84.

### Removed
- The `chocophlan` and `uniref` parameters (now bundle fields) and the
  `withName` entries for `humann` and the HUMAnN DB processes in
  `conf/base.config` (resources are set in the process bodies).

### Changed
- **HUMAnN database caches are keyed by bundle** (`chocophlan/<bundle>/`,
  `uniref/<bundle>/`, `utility_mapping/<bundle>/`) and the DB processes run in
  the bundle's containers rather than the base image; previously cached
  `chocophlan/full/` etc. are not reused (re-download once). ADR-0016
  supersedes ADR-0002.
- ChocoPhlAn, UniRef and utility-mapping databases are fetched only when
  `skip_humann=false`; previously they ran on every run regardless.
- **Reference-database cache (`storeDir`) paths are now keyed by version**
  (#83). Previously every database lived at a fixed name (`metaphlan`,
  `kraken_db`, `card_db`, `card_kma_db`, …), so changing `metaphlan_index`,
  `kraken_db_url` or `card_db_url` silently reused whatever was cached. Layout
  is now `<store_dir>/<db_name>/<version_key>/` (MetaPhlAn: the index; Kraken2
  and CARD: URL basename without archive extension; `card_kma_db` shares the
  CARD key; ChocoPhlAn/UniRef: their params; utility mapping: `full`). The
  KneadData `human_genome`/`mouse_C57BL` paths are unchanged (no version
  parameter selects them). See the `databases.nf` header and the README
  "Reference Database Cache Layout".
- `metaphlan_unknown_viruses_lists`, `metaphlan_unknown_list` and
  `metaphlan_markers` now pass the staged `${metaphlan_db}` to `--db_dir`
  instead of the hardcoded literal `metaphlan`.

### Migration (one-time, before the first run of this version)
Existing stores must be moved into the keyed layout or the databases will be
re-downloaded. Defaults assumed (substitute the key if a different
index/URL was used; skip directories that do not exist). `.command.*` and
`versions.yml` are copied because `storeDir` requires every declared output in
the keyed directory. Verified against a simulated old-layout store with a stub
run (the four databases below were reported as stored/skipped). See the README
"Reference Database Cache Layout" for HUMAnN directories.

```sh
STORE=/projects/seda0001_amc/cmgd/store          # Alpine
STORE=/anvil/projects/x-cis240955/cmgd/store     # Anvil (keep only one STORE= line)

cd "$STORE"
mkdir -p metaphlan.new/mpa_vJan25_CHOCOPhlAnSGB_202503
mv metaphlan metaphlan.new/mpa_vJan25_CHOCOPhlAnSGB_202503/metaphlan
cp -p .command.* versions.yml metaphlan.new/mpa_vJan25_CHOCOPhlAnSGB_202503/
mv metaphlan.new metaphlan

mkdir -p kraken_db.new/k2_pluspf_16_GB_20260226
mv kraken_db kraken_db.new/k2_pluspf_16_GB_20260226/kraken_db
cp -p .command.* kraken_db.new/k2_pluspf_16_GB_20260226/
mv kraken_db.new kraken_db

mkdir -p card_db.new/broadstreet-v4.0.1
mv card_db card_db.new/broadstreet-v4.0.1/card_db
cp -p .command.* card_db.new/broadstreet-v4.0.1/
mv card_db.new card_db

mkdir -p card_kma_db.new/broadstreet-v4.0.1
mv card_kma_db card_kma_db.new/broadstreet-v4.0.1/card_kma_db
cp -p .command.* versions.yml card_kma_db.new/broadstreet-v4.0.1/
mv card_kma_db.new card_kma_db
```

## [2.2.1] - 2026-07-04

### Added
- **`r2` storage profile** publishes outputs to Cloudflare R2
  (`s3://cmgd-raw/cmgd_nextflow/<version>/<sample_id>/…`) and becomes the
  production storage profile (`-profile alpine,r2`), replacing `gcs` per
  nextflow_telemetry ADR 0008. Keys come from `R2_ACCOUNT_ID`,
  `R2_ACCESS_KEY_ID`, `R2_SECRET_ACCESS_KEY`. Needs Nextflow >= 25.04
  (Nextflow 24.04's nf-amazon ignores the R2 endpoint); Alpine now runs
  25.10.8. Output contract unchanged: a test sample reproduced its GCS output
  file-for-file. See ADR-0015 and #79.

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
- **`manifest.json` `software_versions` was incomplete.** It only listed the
  three tools whose `versions.yml` `sample_manifest` staged (fasterq-dump/awscli/
  fastqc from read acquisition, kneaddata/trimmomatic/bowtie2, metaphlan) and
  silently omitted every step added since: `kraken2`, `bracken`, the resistome
  `kma`, the post-decontamination `fastqc`, and `humann` when enabled. Each
  step's `versions.yml` is now collated per sample (via `collectFile`) into the
  single file `build_manifest.py` globs, so `software_versions` reflects every
  tool that actually ran. Skipped optional steps contribute nothing (their
  version channels start empty), and the HUMAnN invocation moved ahead of the
  manifest so its versions are captured when that branch is on.
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
