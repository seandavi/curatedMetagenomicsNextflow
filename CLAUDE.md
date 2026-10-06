# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a Nextflow pipeline for processing human gut metagenomics data. It takes SRA accessions (or local FASTQs), performs host decontamination (KneadData), and produces taxonomic profiles (MetaPhlAn 4.2.2). HUMAnN functional profiling runs as a version-pinned bundle (ADR-0016) with its own MetaPhlAn pass, and is disabled by default (`skip_humann = true`) pending real-sample validation and cost data.

## Commands

All development tasks are driven by `just` recipes (see `justfile`):

```sh
just test-config-local    # Validate Nextflow syntax (local profile)
just test-config-all      # Validate syntax across all profiles
just test-stub            # Fast structural test via -stub-run (~30s)
just test-nf              # Run nf-test suite (unit + smoke tests)
just test-nf-verbose      # Same with detailed output
just test-clean           # Remove test artifacts
just test-bootstrap       # One-time: install nf-test locally
```

To run the pipeline:
```sh
nextflow run main.nf -profile local --metadata_tsv samples.tsv
```

**nf-test does not support running a single test by name in v0.9.0.** Edit individual test files or use `just test-nf-verbose` and grep the output.

## Architecture

### Entry point and module layout

`main.nf` wires channels and orchestrates the workflow. All processes live in `modules/processes/`:

| File | Processes |
|------|-----------|
| `preprocessing.nf` | `fasterq_dump`, `local_fastqc` |
| `profiling.nf` | `kneaddata`, `metaphlan_*`, `sample_to_markers` |
| `rarefaction.nf` | `rarefy_fastq` (seqtk downsampling) |
| `kraken.nf` | `kraken2`, `bracken` (complementary read-based profiling) |
| `resistome.nf` | `resistome_kma` (KMA/CARD read-based AMR profiling, both branches) |
| `qc.nf` | `fastqc` (post-decontamination FastQC) |
| `manifest.nf` | `sample_manifest` (per-sample provenance + read-accounting `manifest.json`) |
| `databases.nf` | Reference database downloads (MetaPhlAn, KneadData, HUMAnN, Kraken2, CARD DBs) |
| `humann.nf` | `metaphlan_for_humann` + `humann` in the selected bundle's containers (disabled by default); wired by the `HUMANN` subworkflow, `modules/subworkflows/humann.nf` |
| `finalize.nf` | `MARK_COMPLETE` sentinel |

DB processes are invoked only from the `DATABASES` subworkflow (`modules/subworkflows/databases.nf`), which gates them per feature (Kraken2 unless `skip_kraken`, CARD/KMA unless `skip_resistome`, HUMAnN DBs only when `!skip_humann`) and emits the DB channels `main.nf` consumes. HUMAnN DBs (bundle-pinned MetaPhlAn index, ChocoPhlAn, UniRef, utility mapping) run in the selected bundle's containers. `--databases_only` (`params.databases_only`) runs just `DATABASES` — no input contract, no per-sample process — to pre-stage `store_dir`.

`sample_manifest` writes one `manifest.json` at each sample's published root via `bin/build_manifest.py` (pure Python, no extra container). It compiles provenance, raw-vs-decontaminated read accounting, rarefaction parameters, and the consolidated per-process `versions.yml`. `MARK_COMPLETE` is gated on it so a sample directory is never marked complete before its manifest exists.

`kraken2`/`bracken` (module `kraken.nf`, gated by `skip_kraken`) add a complementary read-based profile under each branch's `kraken/` subdirectory. They use per-process StaPH-B containers and set all directives (container, resources, `maxForks`) **in the process body** rather than `conf/base.config`, because they are imported under aliases. The Kraken2 DB is loaded into RAM (default mode, not memory-mapped or staged to scratch) for cluster portability; `kraken_maxforks` throttles concurrent reads of the DB off shared storage. See ADR-0006.

`resistome_kma` (module `resistome.nf`, gated by `skip_resistome`) maps the decontaminated reads against a **KMA-indexed CARD** reference (`kma ... -ef`), publishing `card_kma.*.gz` outputs under a `resistome/` subdirectory. CARD is downloaded and its homolog-model FASTA indexed once into `store_dir` (`card_db` → `card_kma_db` in `databases.nf`; the index step runs in the KMA biocontainer). Unlike the former RGI step it runs on **both branches** (imported under `resistome_kma_full`/`resistome_kma_rarefied` aliases), with container/resources set in the process body. The KMA command is not exercised by stub tests; validate it on a real sample. See ADR-0012 (supersedes ADR-0007).

Reference databases are cached under `store_dir` as `<store_dir>/<db_name>/<version_key>/` (the keyed directory holds the `<db_name>` dir plus that task's `.command*`/`versions.yml`; e.g. `metaphlan/<metaphlan_index>/`, `kraken_db/<url-basename-without-extension>/`, `card_db/…`, `card_kma_db/…`; KneadData `human_genome`/`mouse_C57BL` are unversioned). The key convention is documented in the `databases.nf` header; new `storeDir` database processes must follow it, and downstream processes must use the staged input path rather than a hardcoded directory name.

`fastqc` (module `qc.nf`, gated by `skip_fastqc`) runs FastQC on the **decontaminated** reads (`<sample>/fastqc/`); raw-read FastQC already runs in `fasterq_dump`/`local_fastqc`, so this gives a before/after view. It runs in the base image (no new container). A per-sample MultiQC report was considered but dropped to avoid introducing another container (see ADR-0008). Per-sample only.

`HUMANN` subworkflow (full-depth branch only, gated by `!skip_humann`) runs `metaphlan_for_humann` (MetaPhlAn 4.1.1, vJun23 index for the `humann3.9` bundle; independent of the main 4.2.2 pass, whose taxonomy and published paths are untouched) → `humann --taxonomic-profile`. Everything version-dependent — containers, MetaPhlAn version/index, ChocoPhlAn/UniRef/utility-mapping DB names — is one bundle in `conf/humann_bundles.config`, selected by `params.humann_bundle`; an unknown name fails at start (only when HUMAnN is enabled). Outputs publish to `<sample>/humann/<bundle>/` (profile in `metaphlan/`) with HUMAnN-native filenames; container/resources/`maxForks` are in the process body; manifest versions use the distinct keys `metaphlan_humann`/`bowtie2_humann`/`humann` plus `humann_bundle`. HUMAnN DB caches are keyed by bundle. See ADR-0016.

### Architecture Decision Records

Significant, non-obvious decisions are recorded as ADRs in `docs/adr/` (index in `docs/adr/README.md`). Consult them before changing container strategy, the HUMAnN default, the output layout, the manifest, or the Kraken2/resistome steps — and add a new ADR (do not edit accepted ones) when making a comparably significant decision.

### Configuration layers

```
nextflow.config          ← parameters, metadata, reporting
conf/base.config         ← shared resource baselines and retry logic
conf/profiles/           ← composable profiles for compute and storage
  local.config           ← Docker, local scheduler
  google.config          ← Google Batch compute (pair with gcs)
  gcs.config             ← GCS storage (publish_base_dir, workDir, google.project) — legacy
  r2.config              ← Cloudflare R2 storage (publish_base_dir s3://cmgd-raw; keys from R2_* env) — production
  anvil.config           ← SLURM (Anvil HPC)
  alpine.config          ← PBS Pro (CU Alpine HPC)
  unitn.config           ← UNITN HPC
conf/test/               ← overrides for offline/stub testing
nf-test.config           ← nf-test defaults (uses local + stub)
```

Profiles are composable — compute and storage concerns are separated so they can be combined independently:

```sh
-profile google,gcs     # Google Batch compute + GCS output
-profile alpine,r2      # HPC (SLURM) compute + R2 output (production; Nextflow >= 25.04)
-profile alpine,gcs     # HPC (SLURM) compute + GCS output (legacy)
-profile local          # local Docker, output to ./results
```

The default `publish_base_dir` is `${launchDir}/results`. The `r2` profile overrides this to `s3://cmgd-raw` and reads its keys from `R2_ACCOUNT_ID` / `R2_ACCESS_KEY_ID` / `R2_SECRET_ACCESS_KEY` (ADR-0015). The legacy `gcs` profile overrides it to `gs://cmgd-data/results/cMDv<version>` and also sets `google.project`.

### Resource / retry policy (`conf/base.config`)

Memory scales with retries: `effective_memory = baseline × task.attempt`. Exit codes 137–140 (OOM/kill) trigger up to 3 retries. Process labels (`qc`, `profiling`, `db_setup`) map to resource tiers.

### Input modes

- **SRA batch**: `--metadata_tsv samples.tsv` (TSV with `sample_id` + `run_ids` columns)
- **Single sample**: `--sample_id SRS123 --run_ids SRR456`
- **Local FASTQs**: `--local_input true` (skips fasterq_dump)

### Output structure

Outputs follow the pattern: `<publish_base_dir>/<workflow_name>/<pipeline_version>/<sample_id>/<step>/`

### Containers

The base image covers the core toolchain. It is built via `docker/cloudbuild.yaml` (Google Cloud Build, manual trigger), and its tag is pinned in `conf/base.config` under `process.container`. Newer tools run in pinned upstream biocontainers set in their process bodies (ADR-0001). When no usable upstream image exists, the image lives in `docker/<name>/` and a GitHub Actions workflow publishes it to `ghcr.io/seandavi/<name>` (ADR-0017); today that is `docker/humann4a` for the `humann4.0.0a1` bundle.

## Testing

nf-test files live in `tests/`:
- `main.nf.test` — pipeline-level smoke tests
- `main.functions.nf.test` — unit tests for Nextflow functions

Tests use `-stub-run` (no real containers or data needed). The `conf/test/disable-telemetry.config` suppresses outbound calls during CI.
