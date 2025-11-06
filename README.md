# seandavi/curatedmetagenomicsnextflow

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**seandavi/curatedmetagenomicsnextflow** is a bioinformatics pipeline for processing metagenomic sequencing data. The pipeline performs quality control, taxonomic profiling, and optional functional profiling following the curatedMetagenomics workflow.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.

## Pipeline summary

1. Download raw sequencing data from SRA using `fasterq-dump` (or use local FASTQ files)
2. Quality control with `FastQC`
3. Remove contaminating sequences with `KneadData`
4. Taxonomic profiling with `MetaPhlAn`:
   - Generate taxonomic profiles
   - Extract marker information
   - Generate StrainPhlAn marker files
5. Functional profiling with `HUMAnN` (optional):
   - Generate gene family abundance tables
   - Generate pathway abundance and coverage tables
   - Normalize tables (CPM and relative abundance)
   - Split stratified/unstratified tables

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=23.04.0`)

2. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can also use [`Podman`](https://podman.io/))

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run seandavi/curatedmetagenomicsnextflow -profile test,docker --outdir results
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.
   - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.

4. Start running your own analysis!

   ```bash
   nextflow run seandavi/curatedmetagenomicsnextflow --input samplesheet.tsv --outdir results -profile docker
   ```

## Documentation

### Input Samplesheet

You will need to create a samplesheet with information about the samples you would like to analyze before running the pipeline. Use this parameter to specify its location:

```bash
--input '[path to samplesheet file]'
```

The samplesheet should be a tab-separated file with the following columns:

- `sample_id`: Unique sample identifier
- `NCBI_accession`: SRA accession number(s), separated by semicolons for multiple runs

For local FASTQ files (with `--local_input`):

- `sample_id`: Unique sample identifier
- `file_paths`: Path(s) to FASTQ file(s), separated by semicolons for multiple files

Example samplesheet for SRA download:

```
sample_id	NCBI_accession
sample1	SRR1234567
sample2	SRR2345678;SRR2345679
```

Example samplesheet for local files:

```
sample_id	file_paths
sample1	/path/to/sample1_R1.fastq.gz;/path/to/sample1_R2.fastq.gz
sample2	/path/to/sample2_R1.fastq.gz
```

### Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run seandavi/curatedmetagenomicsnextflow --input samplesheet.tsv --outdir results -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Profiles

The pipeline comes with several profiles to suit different execution environments:

- `test`: A minimal test dataset to check pipeline functionality
- `docker`: Use Docker containers
- `singularity`: Use Singularity containers
- `local`: Run locally (requires all software installed)
- `google`: Run on Google Cloud Batch
- `anvil`: Run on AnVIL
- `alpine`: Run on Alpine HPC
- `unitn`: Run on UNITN PBS Pro cluster

You can specify profiles using the `-profile` parameter:

```bash
nextflow run seandavi/curatedmetagenomicsnextflow -profile test,docker
```

Multiple profiles can be specified by separating them with a comma.

### Main arguments

#### `--input`

Path to input samplesheet (TSV format). This replaces the older `--metadata_tsv` parameter.

#### `--outdir`

The output directory where the results will be saved. You must use absolute paths to storage on Cloud infrastructure.

#### `--local_input`

Set to `true` to provide local FASTQ file paths instead of downloading from SRA.

Default: `false`

#### `--skip_humann`

Skip HUMAnN functional profiling step.

Default: `false`

### Reference databases

The pipeline will automatically download and cache reference databases in the location specified by `--store_dir`. These databases include:

- MetaPhlAn database for taxonomic profiling
- ChocoPhlAn pangenome database for HUMAnN
- UniRef protein database for HUMAnN
- Human and mouse reference genomes for KneadData

#### `--store_dir`

Directory to store reference databases.

Default: `'databases'`

#### `--metaphlan_index`

MetaPhlAn database index version.

Default: `'latest'`

#### `--chocophlan`

ChocoPhlAn database version for HUMAnN.

Default: `'full'`

#### `--uniref`

UniRef database version for HUMAnN.

Default: `'uniref90_diamond'`

#### `--organism_database`

Organism reference database for KneadData contamination removal.

Options: `'human_genome'`, `'mouse_C57BL'`

Default: `'human_genome'`

## Output

Results are organized by sample in the output directory:

```
results/
├── sample1/
│   ├── fasterq_dump/          # Raw data download information
│   ├── kneaddata/              # Quality control results
│   ├── metaphlan_lists/        # Taxonomic profiles
│   ├── metaphlan_markers/      # Marker abundance/presence
│   ├── strainphlan_markers/    # Strain-level markers
│   └── humann/                 # Functional profiles (if enabled)
├── sample2/
│   └── ...
└── pipeline_info/              # Pipeline execution information
```

## Credits

seandavi/curatedmetagenomicsnextflow was originally written by Sean Davis.

## Citations

If you use seandavi/curatedmetagenomicsnextflow for your analysis, please cite the following papers:

### Pipeline tools

- [MetaPhlAn](https://github.com/biobakery/MetaPhlAn)

  > Blanco-Míguez A, Beghini F, Cumbo F, et al. Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. Nat Biotechnol. 2023;41(11):1633-1644. doi:10.1038/s41587-023-01688-w

- [HUMAnN](https://github.com/biobakery/humann)

  > Beghini F, McIver LJ, Blanco-Míguez A, et al. Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. eLife. 2021;10:e65088. doi:10.7554/eLife.65088

- [KneadData](https://github.com/biobakery/kneaddata)

  > The KneadData tool is part of the bioBakery suite of tools for metagenomic analysis.

- [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

  > Di Tommaso P, Chatzou M, Floden EW, et al. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017;35(4):316-319. doi:10.1038/nbt.3820

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)
  > Merkel D. Docker: lightweight linux containers for consistent development and deployment. Linux Journal. 2014;2014(239):2.

## Support

For questions or issues, please open an issue on the [GitHub repository](https://github.com/seandavi/curatedmetagenomicsnextflow/issues).
