# Curated Metagenomics NextFlow Pipeline

A NextFlow pipeline for processing metagenomics data, implementing the curatedMetagenomics workflow.

## Overview

This pipeline processes raw sequencing data through multiple steps:
1. FASTQ extraction with `fasterq-dump`
2. Quality control with `KneadData`
3. Taxonomic profiling with `MetaPhlAn`
4. Functional profiling with `HUMAnN` (optional)

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
nextflow run main.nf --metadata_tsv samples.tsv --skip_humann --publish_dir results
```

## Parameters

### General Pipeline Parameters

| Parameter      | Description                            | Default       |
| -------------- | -------------------------------------- | ------------- |
| `metadata_tsv` | Path to TSV file with sample metadata  | `null` |
| `sample_id` | Sample identifier for single-sample mode | `null` |
| `run_ids` | Semicolon-delimited run accessions for single-sample mode | `null` |
| `local_input` | Interpret TSV `file_paths` instead of SRA accessions | `false` |
| `publish_dir`  | Directory to publish results           | `gs://cmgd-data/results/cMDv4` |
| `store_dir`    | Directory to store reference databases | `databases`   |
| `cmgd_version` | Curated Metagenomic Data version       | `4`           |
| `publish_mode` | `publishDir` mode for all published outputs | `copy` |

### Process Control Parameters

| Parameter     | Description                      | Default |
| ------------- | -------------------------------- | ------- |
| `skip_humann` | Skip HUMAnN functional profiling | `true` |

### MetaPhlAn Parameters

| Parameter         | Description            | Default  |
| ----------------- | ---------------------- | -------- |
| `metaphlan_index` | MetaPhlAn index to use | `mpa_vJan25_CHOCOPhlAnSGB_202503` |
| `organism_database` | KneadData reference database | `human_genome` |

### HUMAnN Parameters

| Parameter    | Description                 | Default            |
| ------------ | --------------------------- | ------------------ |
| `chocophlan` | ChocoPhlAn database version | `full`             |
| `uniref`     | UniRef database version     | `uniref90_ec_filtered_diamond` |

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

Results will be organized by sample in the `publish_dir` directory:
```
results/
├── sample1/
│   ├── fasterq_dump/
│   ├── kneaddata/
│   ├── metaphlan_lists/
│   ├── metaphlan_markers/
│   ├── strainphlan_markers/
│   └── humann/
├── sample2/
│   └── ...
```

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
