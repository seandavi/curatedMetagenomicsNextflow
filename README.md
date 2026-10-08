# Curated Metagenomics NextFlow Pipeline

A NextFlow pipeline for processing metagenomics data, implementing the curatedMetagenomics workflow.

## Overview

This pipeline processes raw sequencing data through multiple steps:
1. FASTQ extraction with `fasterq-dump`
2. Quality control with `KneadData` (plus FastQC on the decontaminated reads, optional)
3. Rarefaction to 1 M reads with `seqtk sample` (optional, enabled by default)
4. Taxonomic profiling with `MetaPhlAn` — run on both the full and rarefied reads in parallel
5. Complementary read-based taxonomic profiling with `Kraken2` + `Bracken` (optional, both branches)
6. Antimicrobial-resistance profiling with `KMA` against `CARD` (optional, both branches)
7. Functional profiling with `HUMAnN` (optional, uses full-depth reads)

Each sample also gets a `manifest.json` (provenance + read accounting) and a
`MARK_COMPLETE` sentinel once all enabled branches finish.

## Status and releases

| Line | Where | What it is |
| ---- | ----- | ---------- |
| **2.2.x** | tags `2.2.1`–`2.2.3` on `release/2.2.x` | Current production code. `2.2.3` adds the `r2` storage profile; `2.2.2` points telemetry at the v2 orchestrator. Patch releases only. |
| **2.3.0** (unreleased) | `main` | Next output epoch, for the first full corpus run: MetaPhlAn 4.2.6 with the vJan26 index by default (`mpa4.2.6_vJan26`); HUMAnN only through version-pinned bundles (off by default; no HUMAnN in the base image). Breaking: `--metaphlan_profile` replaces `--metaphlan_index`. Also `--databases_only`, keyed database caches, ENA-first reads. See [`CHANGELOG.md`](CHANGELOG.md). |

The git tag, `manifest.version` in `nextflow.config` and the revision the
orchestrator dispatches move in lockstep. `release/2.2.x` is merged into `main`, so
`main` carries the 2.2.2/2.2.3 changes (telemetry URL, `r2` profile).

## How production runs

Production runs are not launched by hand. The
[nextflow_telemetry](https://github.com/seandavi/nextflow_telemetry) orchestrator
keeps a daemon on each cluster's login node (CU Alpine, Purdue Anvil) that claims
batches of samples and submits one Nextflow driver per batch, roughly:

```bash
nextflow run seandavi/curatedMetagenomicsNextflow -revision <tag> \
  -profile <alpine|anvil>,r2 -c nextflow_override.config \
  --metadata_tsv metadata.tsv --run_name <run> \
  --publish_dir s3://cmgd-raw/<workflow_id>/<version> \
  [-params-file params.json] -with-weblog <telemetry url>
```

- **A registration is one pipeline configuration** (nextflow_telemetry ADR-0010):
  a `workflow_id` + `version` with pinned params, e.g. `cmgd_humann3.9 2.3.0`
  (`skip_humann=false`, `humann_bundle=humann3.9`), `cmgd_humann4a1 2.3.0`
  (pilot collections only) and `cmgd_mpa4.2 2.3.0` (`skip_humann=true`). The
  current production registration is `cmgd_nextflow 2.2.1` running revision `2.2.3`.
- **Results are keyed by registration, not by `manifest.version`**: the
  orchestrator passes `--publish_dir <base>/<workflow_id>/<version>`, so a patch
  revision keeps its registration's prefix and two bundles of one tag never share one.
- **Sample folder names** are the `sample_id` the orchestrator writes into
  `metadata.tsv`: a readset id (`RS.<digest>`, refget seqcol over the run
  accessions; nextflow_telemetry ADR-0007) for registrations created from 2.3.0
  on, the older md5 sample id for `cmgd_nextflow 2.2.1`.
- Telemetry: `api_url` and the weblog URL point at the orchestrator; `run_name`
  is injected by the orchestrator and is `null` in ad-hoc runs (from 2.2.2,
  task-log uploads are then skipped).

## Accessing results

- **Per-sample outputs** are publicly readable (no listing) at
  `https://cmgd-raw.cancerdatasci.org/<workflow_id>/<version>/<sample>/<path>`,
  e.g. `…/cmgd_nextflow/2.2.1/<sample>/manifest.json`. HUMAnN gene families are
  distributed this way, indexed per study in the public releases.
- **Curated tables** (MetaPhlAn, Bracken, resistome, QC, markers, HUMAnN pathways)
  are loaded into a DuckLake and published as versioned, read-only releases at
  `https://cmgd-public.cancerdatasci.org` (DuckDB/Parquet over HTTPS, per-study
  TSV/Parquet downloads). Consumer guide:
  [nextflow_telemetry `docs/data-access.md`](https://github.com/seandavi/nextflow_telemetry/blob/main/docs/data-access.md);
  R client for the cMD team: `cmgdr`.

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

### Stage databases before a batch

Reference databases (MetaPhlAn, KneadData, Kraken2, CARD/KMA, and HUMAnN's
ChocoPhlAn/UniRef/utility mapping when `--skip_humann false`) are downloaded
into `store_dir` on first use, which can take hours. Pre-stage them with
`--databases_only`, which needs no sample inputs, runs only the database
processes (no per-sample step) and exits 0:

```bash
nextflow run main.nf -profile alpine,r2 --databases_only --store_dir /scratch/alpine/$USER/cmgd_db
```

The same skip flags apply, so only the databases a normal run would need are
staged (`--skip_kraken`, `--skip_resistome`, `--skip_humann`). Run it as its
own job (for example an `sbatch` script wrapping the command above) with a
time limit long enough for the downloads, and use the same `store_dir` for the
later batch run: a second invocation finds every database in `store_dir` and
runs no database task.

## Parameters

### General Pipeline Parameters

| Parameter      | Description                            | Default       |
| -------------- | -------------------------------------- | ------------- |
| `metadata_tsv` | Path to TSV file with sample metadata  | `null` |
| `sample_id` | Sample identifier for single-sample mode | `null` |
| `run_ids` | Semicolon-delimited run accessions for single-sample mode | `null` |
| `local_input` | Interpret TSV `file_paths` instead of SRA accessions | `false` |
| `databases_only` | Stage reference databases into `store_dir` and exit (no sample inputs needed) | `false` |
| `publish_base_dir`  | Base directory prefix for published results (`r2` profile: `s3://cmgd-raw`; `gcs` profile: `gs://cmgd-data/results/cMDv4`) | `${launchDir}/results` |
| `publish_dir`  | Full publish root, replacing `<publish_base_dir>/<manifest.name>/<manifest.version>`; the orchestrator sets it to `<base>/<workflow_id>/<version>` | `null` |
| `api_url`      | Orchestrator API base for task-log uploads (paired with the `weblog` URL) | see `nextflow.config` |
| `run_name`     | Orchestrator run name; `null` in ad-hoc runs | `null` |
| `store_dir`    | Directory to store reference databases | `databases`   |
| `cmgd_version` | Curated Metagenomic Data version       | `4`           |
| `publish_mode` | `publishDir` mode for all published outputs | `copy` |

### Reference Database Cache Layout

Reference databases are cached under `store_dir` with Nextflow `storeDir`,
which reuses a task's outputs whenever they already exist. Each
parameter-dependent database therefore lives in a directory keyed by the
version that selects it, `<store_dir>/<db_name>/<version_key>/` (holding the
`<db_name>` directory itself plus that task's `.command*` and `versions.yml`):

| Cache path under `store_dir`                       | Version key                                        |
| -------------------------------------------------- | -------------------------------------------------- |
| `metaphlan/<index>/`                               | the `index` of `metaphlan_profile` (and of a HUMAnN bundle's profile) |
| `kraken_db/<key>/`                                 | `kraken_db_url` basename without archive extension |
| `card_db/<key>/`, `card_kma_db/<key>/`             | `card_db_url` basename without archive extension   |
| `chocophlan/<bundle>/`, `uniref/<bundle>/`, `utility_mapping/<bundle>/` | `humann_bundle` (the bundle pins the DB names) |
| `human_genome/`, `mouse_C57BL/`                    | not versioned (no parameter selects a release)     |

Changing `metaphlan_profile` to one with a different index, `kraken_db_url` or `card_db_url` now creates a new
directory beside the old one rather than silently reusing it. For example, the
2.3.0 default profile installs into `metaphlan/mpa_vJan26_CHOCOPhlAnSGB_202605/`
(about 48 GB to download; stage it with `--databases_only`) next to the 2.2.x
`metaphlan/mpa_vJan25_CHOCOPhlAnSGB_202503/`, and the
default Kraken2 URL `.../k2_pluspf_16_GB_20260226.tar.gz` caches to
`kraken_db/k2_pluspf_16_GB_20260226/` and the default CARD URL
`.../broadstreet-v4.0.1.tar.bz2` to `card_db/broadstreet-v4.0.1/` and
`card_kma_db/broadstreet-v4.0.1/`.

**One-time migration of an existing store.** Stores populated by earlier
versions hold the unkeyed directories (`metaphlan/`, `kraken_db/`, `card_db/`,
`card_kma_db/`) directly in `store_dir`, with a shared `.command.*` and
`versions.yml`. To reuse them instead of re-downloading, move each into its
keyed location once, before the first run of this version. `storeDir` only
hits when every declared output (including `.command*` and, where declared,
`versions.yml`) is present in the keyed directory, so those are copied too.
These commands assume the store was built with the default parameters (if a
different `metaphlan_profile` index / `kraken_db_url` / `card_db_url` was used,
substitute that key); drop the block for any directory that does not exist.
`human_genome/` and `mouse_C57BL/` need no change. The
HUMAnN directories (`chocophlan/`, `uniref/`, `utility_mapping/`) are now keyed
by HUMAnN bundle (default `humann3.9`) rather than by database name, because the
same name can hold different content under different HUMAnN releases. They are
not migrated: if HUMAnN was ever enabled, the databases are downloaded again
into `chocophlan/humann3.9/` etc. (HUMAnN's bundle-pinned MetaPhlAn vJun23 index
is stored under `metaphlan/mpa_vJun23_CHOCOPhlAnSGB_202307/`).

Set `STORE` to the cluster's store (keep only the matching `STORE=` line),
then run the block:

```sh
STORE=/projects/seda0001_amc/cmgd/store          # Alpine
STORE=/anvil/projects/x-cis240955/cmgd/store     # Anvil

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

Verified with a stub run against a simulated old-layout store (the four
databases above were reported as stored and skipped); not run against the real
Alpine/Anvil stores.

### Process Control Parameters

| Parameter        | Description                      | Default |
| ---------------- | -------------------------------- | ------- |
| `skip_humann`    | Skip HUMAnN functional profiling | `true` |
| `skip_rarefied`  | Skip the rarefied profiling branch | `false` |
| `skip_kraken`    | Skip Kraken2 + Bracken read-based profiling | `false` |
| `skip_resistome` | Skip KMA/CARD resistome profiling (both branches) | `false` |
| `skip_fastqc`    | Skip FastQC on the host-decontaminated reads | `false` |

`skip_humann=true` is the default. With `--skip_humann false`, HUMAnN runs on the
full-depth reads in its own version-matched subworkflow (see
[HUMAnN Parameters](#humann-parameters)); enabling it by default is a separate
decision pending real-sample validation and cost data.

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
<sample>/full_data/kraken/
<sample>/rarefied_data/rarefaction/
<sample>/rarefied_data/metaphlan_lists/
<sample>/rarefied_data/metaphlan_markers/
<sample>/rarefied_data/strainphlan_markers/
<sample>/rarefied_data/kraken/
```

Set `--skip_rarefied` to suppress the rarefied branch and restore the original
single-branch layout (without the `full_data/` subdirectory prefix).

### MetaPhlAn Parameters

| Parameter         | Description            | Default  |
| ----------------- | ---------------------- | -------- |
| `metaphlan_profile` | Named MetaPhlAn profile for the main taxonomy pass (see below) | `mpa4.2.6_vJan26` |
| `organism_database` | KneadData reference database | `human_genome` |

A *profile* (defined in [`conf/metaphlan_profiles.config`](conf/metaphlan_profiles.config))
is one pinned unit: container, MetaPhlAn version, index, and the options that
differ between MetaPhlAn releases (`--bowtie2db`/`--bowtie2out` through 4.1.x,
`--db_dir`/`--mapout` from 4.2). An unknown `metaphlan_profile` fails at
start-up with the list of valid names. Profiles:

| Profile | Container | Index |
| ------- | --------- | ----- |
| `mpa4.2.6_vJan26` (default) | `ghcr.io/seandavi/curatedmetagenomics:metaphlan4.2.6` (the base image) | `mpa_vJan26_CHOCOPhlAnSGB_202605` |
| `mpa4.2.2_vJan25` (2.2.x main pass) | `seandavi/curatedmetagenomics:metaphlan4.2.2` (the 2.2.x base image) | `mpa_vJan25_CHOCOPhlAnSGB_202503` |
| `mpa4.1.1_vJun23` | `quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0` | `mpa_vJun23_CHOCOPhlAnSGB_202307` |
| `mpa4.1.1_vOct22` | `quay.io/biocontainers/metaphlan:4.1.1--pyhdfd78af_0` | `mpa_vOct22_CHOCOPhlAnSGB_202403` |

MetaPhlAn 4.2.6 reports itself as `4.2.5` (`--version`; upstream did not bump
it), so that is what `manifest.json` records for the default profile. The
vJan26 SGB-to-GTDB (r226) mapping ships with MetaPhlAn 4.2.5+
(`metaphlan/utils/mpa_vJan26_CHOCOPhlAnSGB_202605_SGB2GTDB_r226.tsv`); GTDB
taxonomy is added in post-processing, not by the pipeline
([ADR-0013](docs/adr/0013-remove-gtdb-conversion.md)).

See [`docs/adr/0018-metaphlan-profiles.md`](docs/adr/0018-metaphlan-profiles.md)
and [`docs/adr/0019-base-image-on-ghcr.md`](docs/adr/0019-base-image-on-ghcr.md).

### Kraken2 / Bracken Parameters

| Parameter             | Description                                                        | Default |
| --------------------- | ------------------------------------------------------------------ | ------- |
| `skip_kraken`         | Disable Kraken2 + Bracken read-based profiling                     | `false` |
| `kraken_db_url`       | Prebuilt Kraken2 index tarball (bundles Bracken kmer distributions) | PlusPF 16 GB cap |
| `kraken_confidence`   | Kraken2 `--confidence` threshold (0.0–1.0)                         | `0.0`   |
| `kraken_maxforks`     | Max concurrent Kraken2 tasks (throttles shared-storage DB reads)   | `4`     |
| `bracken_read_length` | Bracken read length (must have a matching kmer distribution)       | `100`   |

When `skip_kraken=false` (the default), the host-decontaminated reads are
classified with [Kraken2](https://ccb.jhu.edu/software/kraken2/) and re-estimated
with [Bracken](https://ccb.jhu.edu/software/bracken/) at species and genus level,
complementing MetaPhlAn's marker-based view with a whole-community, read-count
profile. Outputs (`kraken2.report.txt.gz`, `bracken.species.txt.gz`,
`bracken.genus.txt.gz`, and the Bracken-adjusted reports) are published under a
`kraken/` subdirectory in each branch.

The default database is the **16 GB capped PlusPF** index (RefSeq bacteria/
archaea/viral + protozoa + fungi + human decoy) — broad enough for a multi-body-
site cohort, including the skin/airway mycobiome. Prebuilt PlusPF indexes come
in 8 GB, 16 GB, or full sizes (no 32 GB); point `kraken_db_url` at the full
tarball for better rare-taxon sensitivity at higher RAM/I/O cost. The database is
downloaded once into `store_dir` and **loaded into RAM** at run time (not
memory-mapped or staged to node-local scratch) for cluster portability;
`kraken_maxforks` limits how many tasks read it off shared storage at once. See
[`docs/adr/0006-kraken2-bracken-complementary-profiler.md`](docs/adr/0006-kraken2-bracken-complementary-profiler.md).

### Resistome (KMA/CARD) Parameters

| Parameter        | Description                                          | Default |
| ---------------- | ---------------------------------------------------- | ------- |
| `skip_resistome` | Disable KMA/CARD resistome profiling                 | `false` |
| `card_db_url`    | CARD "broadstreet" release tarball URL               | `https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2` |

When `skip_resistome=false` (the default), the host-decontaminated reads are
profiled for antimicrobial-resistance genes with [KMA](https://bitbucket.org/genomicepidemiology/kma)
against a KMA-indexed [CARD](https://card.mcmaster.ca/) reference, emitting
per-sample template hits with depth/coverage statistics (`card_kma.res.gz`,
`card_kma.mapstat.gz`, plus `card_kma.{aln,fsa,frag}.gz`) under a `resistome/`
subdirectory. The CARD broadstreet release is downloaded once into `store_dir`,
and its homolog-model nucleotide FASTA is indexed with `kma index` (also cached
in `store_dir`) for reuse across runs.

Unlike the previous RGI/CARD step (full branch only), KMA is cheap enough to run
on **both the full and rarefied branches**, matching MetaPhlAn and Kraken. ARGs
are sparse at the rarefied 1M-read depth, so treat the rarefied resistome as
low-sensitivity. See
[`docs/adr/0012-resistome-kma-card.md`](docs/adr/0012-resistome-kma-card.md).

### QC Parameters

| Parameter     | Description                                          | Default |
| ------------- | ---------------------------------------------------- | ------- |
| `skip_fastqc` | Disable FastQC on the host-decontaminated reads      | `false` |

Raw-read FastQC already runs in both input modes (inside `fasterq_dump` /
`local_fastqc`). When `skip_fastqc=false` (the default), the pipeline also runs
FastQC on the **host-decontaminated** reads, published under `<sample>/fastqc/`,
giving a before/after view of what QC did to the reads. It runs in the base
image, so no additional container is required. Per-sample only. See
[`docs/adr/0008-per-sample-qc-reporting.md`](docs/adr/0008-per-sample-qc-reporting.md).

### HUMAnN Parameters

| Parameter         | Description                                              | Default     |
| ----------------- | -------------------------------------------------------- | ----------- |
| `humann_bundle`   | Named, version-pinned HUMAnN bundle (see below)          | `humann3.9` |
| `humann_maxforks` | Max concurrent `humann` tasks (shared-storage DB reads)  | `4`         |

A *bundle* (defined in [`conf/humann_bundles.config`](conf/humann_bundles.config))
pins the HUMAnN container, the *name* of a MetaPhlAn profile (see above), and
the ChocoPhlAn/UniRef/utility-mapping database names. There are no per-tool
version parameters; an unknown `humann_bundle` (or a bundle naming an unknown
profile) fails at start-up with the list of valid names (checked only when
`--skip_humann false`). Bundles:

| Bundle | HUMAnN image | MetaPhlAn profile |
| ------ | ------------ | ----------------- |
| `humann3.9` (default) | `quay.io/biocontainers/humann:3.9--py312hdfd78af_0` | `mpa4.1.1_vJun23` (HUMAnN 3.9 rejects any other index) |
| `humann4.0.0a1` | `ghcr.io/seandavi/humann:4.0.0a1`, built from [`docker/humann4a`](docker/humann4a/Dockerfile) ([ADR-0017](docs/adr/0017-self-built-images-on-ghcr.md)) | `mpa4.1.1_vOct22` |

HUMAnN 4.0.0a1 is an alpha. Its native output names differ (`out_2_genefamilies`,
`out_3_reactions`, `out_4_pathabundance`, `out_5_pathcoverage`) and it reports
unmapped reads as `READS_UNMAPPED`.

HUMAnN runs on the full-depth branch only. When the bundle's MetaPhlAn profile
differs from `metaphlan_profile` (both current bundles), HUMAnN runs its own
MetaPhlAn pass, so the main-pass taxonomy (MetaPhlAn 4.2.6 / vJan26 by default) is unaffected. When
they are equal (a future HUMAnN release that accepts the main profile), HUMAnN
reuses the main pass's full-depth profile and no second MetaPhlAn or index
install runs. Outputs are published
under `<sample>/humann/<humann_bundle>/` using HUMAnN's native filenames, with
the profile that drove the stratification (a copy of the main profile, on reuse)
in `<sample>/humann/<humann_bundle>/metaphlan/`.
Functional profiles are stratified by that bundle's taxonomy, not by the
`metaphlan_lists`/`metaphlan_markers` profiles. See
[`docs/adr/0016-humann-bundles.md`](docs/adr/0016-humann-bundles.md) and
[`docs/adr/0018-metaphlan-profiles.md`](docs/adr/0018-metaphlan-profiles.md).

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

Results are organized by sample under the publish root: `publish_dir` when set
(production: `<base>/<workflow_id>/<version>`), otherwise
`<publish_base_dir>/<manifest.name>/<manifest.version>`, as in the trees below.

### Dual-branch layout (default, `skip_rarefied=false`)

```
<publish_base_dir>/
├── cmgd_nextflow/
│   ├── 2.3.0/
│   │   ├── sample1/
│   │   │   ├── manifest.json     (provenance + read accounting)
│   │   │   ├── MARK_COMPLETE
│   │   │   ├── kneaddata/
│   │   │   ├── fastqc/           (decontaminated-read FastQC; only when --skip_fastqc false)
│   │   │   ├── full_data/
│   │   │   │   ├── metaphlan_lists/
│   │   │   │   ├── metaphlan_markers/
│   │   │   │   ├── strainphlan_markers/
│   │   │   │   ├── kraken/       (only when --skip_kraken false)
│   │   │   │   └── resistome/    (only when --skip_resistome false)
│   │   │   ├── humann/<bundle>/  (only when --skip_humann false; full-depth reads; profile in metaphlan/)
│   │   │   └── rarefied_data/
│   │   │       ├── rarefaction/
│   │   │       ├── metaphlan_lists/
│   │   │       ├── metaphlan_markers/
│   │   │       ├── strainphlan_markers/
│   │   │       ├── kraken/       (only when --skip_kraken false)
│   │   │       └── resistome/    (only when --skip_resistome false)
│   │   ├── sample2/
│   │   │   └── ...
```

### Single-branch layout (`--skip_rarefied`, backward-compatible)

```
<publish_base_dir>/
├── cmgd_nextflow/
│   ├── 2.3.0/
│   │   ├── sample1/
│   │   │   ├── manifest.json   (provenance + read accounting)
│   │   │   ├── MARK_COMPLETE
│   │   │   ├── kneaddata/
│   │   │   ├── fastqc/         (decontaminated-read FastQC; only when --skip_fastqc false)
│   │   │   ├── metaphlan_lists/
│   │   │   ├── metaphlan_markers/
│   │   │   ├── strainphlan_markers/
│   │   │   ├── kraken/         (only when --skip_kraken false)
│   │   │   ├── resistome/      (only when --skip_resistome false)
│   │   │   └── humann/<bundle>/  (only when --skip_humann false; profile in metaphlan/)
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
- `alpine`: CU Boulder Alpine (SLURM) — production
- `anvil`: Purdue Anvil (SLURM, ACCESS allocation) — production
- `google`: For execution on Google Cloud Batch
- `unitn`: For execution on UNITN PBS Pro

Storage profiles (compose with a compute profile; they only change where outputs are published):
- `r2`: Publish to Cloudflare R2 `s3://cmgd-raw` (production; needs Nextflow ≥ 25.04 and `R2_ACCOUNT_ID`, `R2_ACCESS_KEY_ID`, `R2_SECRET_ACCESS_KEY` in the environment; see ADR-0015)
- `gcs`: Publish to GCS `gs://cmgd-data/results/cMDv4` (legacy; no new production writes)

Production composes a cluster with `r2`, e.g. `-profile alpine,r2` or
`-profile anvil,r2`. Example:
```bash
nextflow run main.nf -profile alpine,r2 --metadata_tsv samples.tsv
```

## Resource And Retry Policy

The baseline CPU, memory and time requests are defined per process in
`conf/base.config` (and, for processes imported under aliases, in the process
body). From 2.3.0 they are sized from real 2.2.x traces
([nextflow_telemetry `docs/research/resource-tuning-2.3.0.md`](https://github.com/seandavi/nextflow_telemetry/blob/main/docs/research/resource-tuning-2.3.0.md)):
memory is set close to each step's measured peak plus headroom, and
single-threaded steps get 1-2 cpus. MetaPhlAn gets 55 GiB, an estimate for the
vJan26 index that the rehearsal should re-measure. On Anvil, the full-branch
MetaPhlAn also gets 30 cpus, which its memory already pays for there. For
retries, memory and time are scaled linearly by retry attempt:

```text
effective_memory = baseline_memory * task.attempt
effective_time   = baseline_time * task.attempt
```

Exit codes 137-140 (OOM or scheduler kill) are retried up to 4 times with that
growth. When the OOM killer takes bowtie2 inside MetaPhlAn, MetaPhlAn itself
exits 1. `bin/bowtie2_oom_to_137` wraps the MetaPhlAn steps that run bowtie2 and
turns that case into 137, so the OOM retry applies (#107).

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

> **Before a full-scale run**, stub tests do not pull containers or download
> databases — work through
> [`docs/verify-before-full-run.md`](docs/verify-before-full-run.md) to confirm
> the pinned container tags and database URLs resolve and the real tool
> invocations behave.

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
- Nextflow ≥ 25.04 for the `r2` profile (production runs 25.10.8), and < 26.04
  until the pipeline passes Nextflow's strict syntax parser
- Java 17+ (Anvil uses a user-space Temurin 21 JDK because its modules stop at Java 11)
- Container support (Singularity/Apptainer on the HPC clusters, Docker locally)
