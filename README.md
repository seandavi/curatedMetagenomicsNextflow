# Curated Metagenomics NextFlow Pipeline

A NextFlow pipeline for processing metagenomics data, implementing the curatedMetagenomics workflow.

## Overview

This pipeline processes raw sequencing data through multiple steps:
1. FASTQ extraction with `fasterq-dump`
2. Quality control with `KneadData`
3. Rarefaction to 1 M reads with `seqtk sample` (optional, enabled by default)
4. Taxonomic profiling with `MetaPhlAn` — run on both the full and rarefied reads in parallel
5. Functional profiling with `HUMAnN` (optional, uses full-depth reads)

## Repository Structure

The pipeline is organized so that workflow orchestration and operational policy are easy to find:

- `main.nf`
  - top-level workflow orchestration and sample-channel wiring
- `modules/processes/`
  - grouped process implementations by functional area
- `nextflow.config`
  - user-facing defaults, metadata, reporting, and config includes
- `conf/base.config`
  - shared retry, container, and per-process resource policy
- `conf/profiles/*.config`
  - site- and executor-specific profile settings
- `conf/test/disable-telemetry.config`
  - offline local smoke-test overrides for `-stub-run`

## Usage

Basic usage:

```bash
nextflow run main.nf --metadata_tsv samples.tsv
```

With specific parameters:

```bash
nextflow run main.nf --metadata_tsv samples.tsv --skip_humann --publish_base_dir results
```

## Parameters

### General Pipeline Parameters

| Parameter      | Description                            | Default       |
| -------------- | -------------------------------------- | ------------- |
| `metadata_tsv` | Path to TSV file with sample metadata  | `null` |
| `sample_id` | Sample identifier for single-sample mode | `null` |
| `run_ids` | Semicolon-delimited run accessions for single-sample mode | `null` |
| `local_input` | Interpret TSV `file_paths` instead of SRA accessions | `false` |
| `publish_base_dir`  | Base directory prefix for published results | `gs://cmgd-data/results/cMDv4` |
| `publish_dir`  | Optional full publish root override after workflow-name/version expansion | `null` |
| `store_dir`    | Directory to store reference databases | `databases`   |
| `cmgd_version` | Curated Metagenomic Data version       | `4`           |
| `publish_mode` | `publishDir` mode for all published outputs | `copy` |

### Process Control Parameters

| Parameter        | Description                      | Default |
| ---------------- | -------------------------------- | ------- |
| `skip_humann`    | Skip HUMAnN functional profiling | `true` |
| `skip_rarefied`  | Skip the rarefied profiling branch | `false` |
| `skip_gtdb`      | Skip MetaPhlAn-to-GTDB profile conversion | `false` |

`skip_humann=true` is the current supported default. The `skip_humann=false`
path is kept in the pipeline for future use, but it is not expected to work
correctly at present because the active MetaPhlAn database/version combination
is not aligned with the HUMAnN branch.

### Rarefaction Parameters

| Parameter      | Description                                       | Default   |
| -------------- | ------------------------------------------------- | --------- |
| `skip_rarefied`| Disable the rarefied branch (restores legacy layout) | `false` |
| `rarefy_reads` | Target read depth for rarefaction                 | `1000000` |
| `rarefy_seed`  | Random seed for reproducible rarefaction          | `42`      |

When `skip_rarefied=false` (the default) the pipeline runs a rarefied profiling
branch in parallel with the full-depth branch.  Both sets of outputs appear
under the sample directory, distinguished by a branch subdirectory:

```
<sample>/full_data/metaphlan_lists/
<sample>/full_data/metaphlan_markers/
<sample>/full_data/strainphlan_markers/
<sample>/full_data/gtdb/
<sample>/rarefied_data/rarefaction/
<sample>/rarefied_data/metaphlan_lists/
<sample>/rarefied_data/metaphlan_markers/
<sample>/rarefied_data/strainphlan_markers/
<sample>/rarefied_data/gtdb/
```

Set `--skip_rarefied` to suppress the rarefied branch and restore the original
single-branch layout (without the `full_data/` subdirectory prefix).

### MetaPhlAn Parameters

| Parameter         | Description            | Default  |
| ----------------- | ---------------------- | -------- |
| `metaphlan_index` | MetaPhlAn index to use | `mpa_vJan25_CHOCOPhlAnSGB_202503` |
| `organism_database` | KneadData reference database | `human_genome` |

### GTDB Parameters

| Parameter      | Description                                          | Default |
| -------------- | ---------------------------------------------------- | ------- |
| `skip_gtdb`    | Disable MetaPhlAn-to-GTDB profile conversion         | `false` |
| `sgb2gtdb_url` | URL of the SGB-to-GTDB assignment table (must match `metaphlan_index`) | MetaPhlAn `mpa_vJan25_CHOCOPhlAnSGB_202503` table |

When `skip_gtdb=false` (the default), each MetaPhlAn relative-abundance profile
is translated into a [GTDB-taxonomy](https://gtdb.ecogenomic.org/) profile using
MetaPhlAn's official SGB-to-GTDB assignment table. The table is downloaded once
into `store_dir` (process `sgb_to_gtdb_db`) and reused across runs. Conversion
is performed by the vendored `bin/sgb_to_gtdb_profile.py` helper, which
substitutes 1:1 mappings directly, bins (sums) n:1 mappings into the shared GTDB
taxon, and aggregates abundances up the GTDB lineage. Output is published as
`gtdb_profile.tsv.gz` under a `gtdb/` subdirectory in each branch.

`sgb2gtdb_url` must be kept in lockstep with `metaphlan_index`: if the MetaPhlAn
index changes, point `sgb2gtdb_url` at the assignment table that ships with the
new index.

### HUMAnN Parameters

| Parameter    | Description                 | Default            |
| ------------ | --------------------------- | ------------------ |
| `chocophlan` | ChocoPhlAn database version | `full`             |
| `uniref`     | UniRef database version     | `uniref90_ec_filtered_diamond` |

Important: enabling the HUMAnN branch currently requires coordinated
MetaPhlAn/HUMAnN version alignment and validation. Treat the existing HUMAnN
path as preserved-but-dormant until that compatibility work is done.

## Input Format

The `metadata_tsv` file should be a tab-separated values file with at least the following columns:
- `sample_id`: Unique sample identifier
- `NCBI_accession`: SRA accession number(s), separated by semicolons for multiple files

If `--local_input true` is used, the TSV should provide:
- `sample_id`
- `file_paths`: Semicolon-delimited local FASTQ paths

Example:
```
sample_id    NCBI_accession
sample1      SRR1234567
sample2      SRR2345678;SRR2345679
```

Local-input example:
```
sample_id    file_paths
sample1      /data/sample1_R1.fastq.gz;/data/sample1_R2.fastq.gz
```

## Output

Results will be organized by sample in the `publish_dir` directory.

### Dual-branch layout (default, `skip_rarefied=false`)

```
<publish_base_dir>/
├── cmgd_nextflow/
│   ├── 1.6.0/
│   │   ├── sample1/
│   │   │   ├── manifest.json     (provenance + read accounting)
│   │   │   ├── MARK_COMPLETE
│   │   │   ├── kneaddata/
│   │   │   ├── full_data/
│   │   │   │   ├── metaphlan_lists/
│   │   │   │   ├── metaphlan_markers/
│   │   │   │   ├── strainphlan_markers/
│   │   │   │   └── gtdb/         (only when --skip_gtdb false)
│   │   │   └── rarefied_data/
│   │   │       ├── rarefaction/
│   │   │       ├── metaphlan_lists/
│   │   │       ├── metaphlan_markers/
│   │   │       ├── strainphlan_markers/
│   │   │       └── gtdb/         (only when --skip_gtdb false)
│   │   ├── sample2/
│   │   │   └── ...
```

### Single-branch layout (`--skip_rarefied`, backward-compatible)

```
<publish_base_dir>/
├── cmgd_nextflow/
│   ├── 1.6.0/
│   │   ├── sample1/
│   │   │   ├── manifest.json   (provenance + read accounting)
│   │   │   ├── MARK_COMPLETE
│   │   │   ├── kneaddata/
│   │   │   ├── metaphlan_lists/
│   │   │   ├── metaphlan_markers/
│   │   │   ├── strainphlan_markers/
│   │   │   ├── gtdb/           (only when --skip_gtdb false)
│   │   │   └── humann/         (only when --skip_humann false)
│   │   ├── sample2/
│   │   │   └── ...
```

### Per-sample manifest

Every sample gets a single `manifest.json` at its published root that compiles:

- **provenance** — pipeline/Nextflow versions, container image, command line,
  run name, git commit, key parameters, and input accessions/paths;
- **read accounting** — raw and host-decontaminated read counts, base counts,
  read-length statistics (min/median/max/mean) and GC%, plus the surviving
  read/base fractions;
- **rarefaction** parameters (when the rarefied branch is active);
- **software_versions** — the per-process `versions.yml` files consolidated
  into one tool→version map.

Read statistics are computed by `bin/build_manifest.py` (pure Python, single
streaming pass, no extra container). `MARK_COMPLETE` is gated on the manifest,
so a sample directory is never marked complete before `manifest.json` exists.
See [`docs/adr/0005-per-sample-manifest.md`](docs/adr/0005-per-sample-manifest.md).

The sample-level directory name is the normalized `meta.sample` value used for
task tags. That comes from `--sample_id` in single-sample mode or the
`sample_id` column in TSV-driven runs.

## Profiles

The pipeline comes with several execution profiles:
- `local`: For local execution
- `google`: For execution on Google Cloud Batch
- `anvil`: For execution on AnVIL
- `alpine`: For execution on Alpine HPC
- `unitn`: For execution on UNITN PBS Pro

Example:
```bash
nextflow run main.nf -profile google --metadata_tsv samples.tsv
```

## Resource And Retry Policy

The baseline CPU and memory requests are defined per process in `conf/base.config`.
Attempt 1 preserves the original resource baselines from the pipeline. For retries,
memory is scaled linearly by retry attempt:

```text
effective_memory = baseline_memory * task.attempt
```

This allows failed tasks to request more memory on later attempts without changing
the initial scheduling footprint.

Process labels in `main.nf` and `modules/processes/*.nf` are semantic rather than
prescriptive. They are intended to make the pipeline easier to read and to support
future policy tuning without hiding the current per-process baselines.

## Local Structural Testing

To validate workflow structure without a cluster or container runtime, use
Nextflow stubs with the offline test override config:

```bash
NXF_DISABLE_CHECK_LATEST=true \
nextflow run . \
  -profile local \
  -stub-run \
  -c conf/test/disable-telemetry.config \
  --sample_id TEST_SAMPLE \
  --run_ids SRR000001
```

That test path disables telemetry, reports, trace output, containers, and cloud
publishing so the full DAG can be validated in restricted local environments.

## Development Test Harness

The recommended developer test entrypoint is `just`, with Nextflow-native tests
implemented in `nf-test`.

Why this split:

- `just` gives the repository one obvious command surface for developers and CI
- `nextflow config` checks catch profile and config composition regressions
- `nextflow ... -stub-run` provides a fast whole-DAG structural smoke test
- `nf-test` is purpose-built for testing Nextflow pipelines, processes, and functions

### `just` commands

This repository currently keeps all developer recipes in a single `test` group:

```bash
just --list --unsorted
```

The most important recipes are:

- `just test-bootstrap`
- `just test-config-local`
- `just test-config-all`
- `just test-stub`
- `just test-nf`
- `just test-clean`

### `nf-test` layout

The repo includes a minimal `nf-test` setup intended for structural and helper
function validation rather than full end-to-end execution:

- `nf-test.config`
  - default nf-test settings
- `tests/nextflow.config`
  - test-specific Nextflow config that reuses the offline local overrides
- `tests/main.nf.test`
  - top-level pipeline smoke tests
- `tests/main.functions.nf.test`
  - helper-function unit tests

`nf-test` is not bundled with Nextflow and it is not a Nextflow plugin. This
repository treats it as a repo-local developer tool and installs it into
`tools/bin/nf-test`.

Bootstrap it with:

```bash
just test-bootstrap
```

After that, use the `just` recipes rather than relying on a separately managed
global `nf-test` binary.

### Recommended local workflow

For routine development, use this order:

1. `just test-bootstrap`
2. `just test-config-local`
3. `just test-stub`
4. `just test-nf`

That sequence is intentionally lightweight. It gives fast feedback on config
resolution, workflow structure, and helper-function behavior without requiring a
cluster or a complete real-data e2e suite.

## Invariants

The following are treated as compatibility constraints during refactoring:

- published subdirectory names must remain unchanged
- published filenames must remain unchanged
- process output filenames must remain unchanged
- sample completion semantics for `MARK_COMPLETE` must remain unchanged

## Dependencies

This pipeline requires:
- Nextflow 22.10.0 or later
- Container support (Docker, Singularity, etc.)
