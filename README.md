# Curated Metagenomics NextFlow Pipeline

A NextFlow pipeline for processing metagenomics data, implementing the curatedMetagenomics workflow.

## Overview

This pipeline processes raw sequencing data through multiple steps:
1. FASTQ extraction with `fasterq-dump`
2. Quality control with `KneadData`
3. Taxonomic profiling with `MetaPhlAn`
4. Functional profiling with `HUMAnN` (optional)

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
| `metadata_tsv` | Path to TSV file with sample metadata  | `samples.tsv` |
| `publish_dir`  | Directory to publish results           | `results`     |
| `store_dir`    | Directory to store reference databases | `databases`   |
| `cmgd_version` | Curated Metagenomic Data version       | `4`           |

### Process Control Parameters

| Parameter     | Description                      | Default |
| ------------- | -------------------------------- | ------- |
| `skip_humann` | Skip HUMAnN functional profiling | `false` |

### MetaPhlAn Parameters

| Parameter         | Description            | Default  |
| ----------------- | ---------------------- | -------- |
| `metaphlan_index` | MetaPhlAn index to use | `latest` |

### HUMAnN Parameters

| Parameter    | Description                 | Default            |
| ------------ | --------------------------- | ------------------ |
| `chocophlan` | ChocoPhlAn database version | `full`             |
| `uniref`     | UniRef database version     | `uniref90_diamond` |

## Input Format

The `metadata_tsv` file should be a tab-separated values file with at least the following columns:
- `sample_id`: Unique sample identifier
- `NCBI_accession`: SRA accession number(s), separated by semicolons for multiple files

Example:
```
sample_id    NCBI_accession
sample1      SRR1234567
sample2      SRR2345678;SRR2345679
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

## Dependencies

This pipeline requires:
- Nextflow 22.10.0 or later
- Container support (Docker, Singularity, etc.)
