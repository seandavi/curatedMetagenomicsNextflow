# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a Nextflow pipeline for processing human gut metagenomics data. It takes SRA accessions (or local FASTQs), performs host decontamination (KneadData), and produces taxonomic profiles (MetaPhlAn 4.2.2). HUMAnN functional profiling is present but disabled by default (`skip_humann = true`) due to unresolved version compatibility issues.

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
| `databases.nf` | Reference database downloads (MetaPhlAn, KneadData, HUMAnN DBs) |
| `humann.nf` | HUMAnN gene/pathway abundance (disabled by default) |
| `finalize.nf` | `MARK_COMPLETE` sentinel |

### Configuration layers

```
nextflow.config          ← parameters, metadata, reporting
conf/base.config         ← shared resource baselines and retry logic
conf/profiles/           ← composable profiles for compute and storage
  local.config           ← Docker, local scheduler
  google.config          ← Google Batch compute (pair with gcs)
  gcs.config             ← GCS storage (publish_base_dir, workDir, google.project)
  anvil.config           ← SLURM (Anvil HPC)
  alpine.config          ← PBS Pro (CU Alpine HPC)
  unitn.config           ← UNITN HPC
conf/test/               ← overrides for offline/stub testing
nf-test.config           ← nf-test defaults (uses local + stub)
```

Profiles are composable — compute and storage concerns are separated so they can be combined independently:

```sh
-profile google,gcs     # Google Batch compute + GCS output
-profile alpine,gcs     # HPC (PBS Pro) compute + GCS output
-profile local          # local Docker, output to ./results
```

The default `publish_base_dir` is `${launchDir}/results`. The `gcs` profile overrides this to `gs://cmgd-data/results/cMDv<version>` and also sets `workDir` and `google.project`.

### Resource / retry policy (`conf/base.config`)

Memory scales with retries: `effective_memory = baseline × task.attempt`. Exit codes 137–140 (OOM/kill) trigger up to 3 retries. Process labels (`qc`, `profiling`, `db_setup`) map to resource tiers.

### Input modes

- **SRA batch**: `--metadata_tsv samples.tsv` (TSV with `sample_id` + `run_ids` columns)
- **Single sample**: `--sample_id SRS123 --run_ids SRR456`
- **Local FASTQs**: `--local_input true` (skips fasterq_dump)

### Output structure

Outputs follow the pattern: `<publish_base_dir>/<workflow_name>/<pipeline_version>/<sample_id>/<step>/`

### Containers

A single Docker image covers all tools. It is built via `docker/cloudbuild.yaml` (Google Cloud Build, manual trigger). The image tag is pinned in `nextflow.config` under `process.container`.

## Testing

nf-test files live in `tests/`:
- `main.nf.test` — pipeline-level smoke tests
- `main.functions.nf.test` — unit tests for Nextflow functions

Tests use `-stub-run` (no real containers or data needed). The `conf/test/disable-telemetry.config` suppresses outbound calls during CI.
