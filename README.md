# README

## Overview

This repo is a fork of [seandavi/curatedMetagenomicsNextflow](https://github.com/seandavi/curatedMetagenomicsNextflow), adapted to work with the University of Trento's HPC cluster and incorporate additional MetaPhlAn outputs.

### Changes

 - Update Dockerfile to include latest versions of Bowtie2 and SAMtools (required for new processes)
 - Add viral profiling, unknown estimation, and consensus-marker processes
 - Enable `.tsv` sample input
 - Establish filepaths specific to UniTn cluster
 - Add `unitn` profile
 - Create submit script `submit_unitn.sh` for UniTn cluster
 - Add documentation for running on the UniTn cluster and handling known errors

## Pipeline

### Original pipeline

- fasterq-dump + FastQC
- setup databases
  - MetaPhlAn
  - ChocoPhlAn
  - UniRef
  - KneadData
    - human genome
    - ribosomal RNA
- MetaPhlAn bug list
- MetaPhlAn markers
- HUMAnN (including various aggregations)

### Added outputs

 - MetaPhlAn viruses list
 - MetaPhlAn unknown estimation
 - MetaPhlAn consensus-marker files

## Inputs and Outputs

### Input table

The `metadata_tsv` file must be:

- tab-separated
- must contain columns
  - `sample_id`
  - `NCBI_accession`, a semicolon-separated list of SRRs
- Can be a file or a web url

If using a Google Bucket, the name bucket must not have underscores.

### Sample IDs

Sample IDs of any sort can be used, but it is recommended to use UUIDs for pipeline submission and maintain a map between them and more readable sample identifiers.

### Output files grouped by tool:

**fasterq-dump + FastQC**

- `fastq_line_count.txt`
- `sampleinfo.txt`
- `fastqc_data.txt`

**KneadData**

- `kneaddata_fastq_linecounts.txt`
- `out_kneaddata.log`

**MetaPhlAn**

- `metaphlan_bugs_list.tsv.gz`
- `metaphlan_viruses_list.tsv.gz`*
- `metaphlan_unknown_list.tsv.gz`*
- `marker_abundance.tsv.gz`
- `marker_presence.tsv.gz`
- `sam.json.bz2`*

=======
\* = new with this fork

**HUMAnN**

- `out_genefamilies.tsv.gz`
- `out_genefamilies_cpm.tsv.gz`
- `out_genefamilies_relab.tsv.gz`
- `out_genefamilies_stratified.tsv.gz`
- `out_genefamilies_unstratified.tsv.gz`
- `out_genefamilies_cpm_stratified.tsv.gz`
- `out_genefamilies_relab_stratified.tsv.gz`
- `out_genefamilies_cpm_unstratified.tsv.gz`
- `out_genefamilies_relab_unstratified.tsv.gz`
- `out_pathabundance.tsv.gz`
- `out_pathabundance_cpm.tsv.gz`
- `out_pathabundance_relab.tsv.gz`
- `out_pathabundance_stratified.tsv.gz`
- `out_pathabundance_unstratified.tsv.gz`
- `out_pathabundance_cpm_stratified.tsv.gz`
- `out_pathabundance_relab_stratified.tsv.gz`
- `out_pathabundance_cpm_unstratified.tsv.gz`
- `out_pathabundance_relab_unstratified.tsv.gz`
- `out_pathcoverage_unstratified.tsv.gz`
- `out_pathcoverage_stratified.tsv.gz`
- `out_pathcoverage.tsv.gz`

## Install

```sh
export NXF_MODE=google
curl https://get.nextflow.io | bash
```

## Google Setup

You will need to be able to access google cloud storage as well as to
run the Google Cloud Pipeline API. This requires credentials to do so.
You can either use Google Default Application Credentials or a key file.
The latter is the recommended approach. If you need a keyfile, contact
the person who owns the Google Project you'll be using.

### Keyfile setup

Once you get a keyfile (which is a json file), run the following:

```sh
# just examples:
export SVC_ACCOUNT='nextflow-service-account@curatedmetagenomicdata.iam.gserviceaccount.com' #example name
export GOOGLE_APPLICATION_CREDENTIALS=/data/curatedmetagenomicdata-f7fc1489b036.json
export GCP_PROJECT=curatedmetagenomicdata
gcloud auth activate-service-account \
   $SVC_ACCOUNT \
   --key-file=$GOOGLE_APPLICATION_CREDENTIALS \
   --project=$GCP_PROJECT
```

## Execution

This assumes that you are running on Google, that credentials are set up correctly,
and that you have a Google Storage Bucket already created. Note that bucket names
must NOT contain the `_` or other special characters.

```sh
# No '_' or other non-url-safe characters in bucket names
export GOOGLE_BUCKET_NAME='your-bucket-name'

# if bucket does not exist:
gsutil mb gs://$GOOGLE_BUCKET_NAME
```

You can now run test data. This will take a few hours the first time, so run on a system that will remain on
during that time (laptops are not a good choice if you are going to close it and go home, for example).

```sh
./nextflow run seandavi/curatedMetagenomicsNextflow \
  -profile google \
  -work-dir gs://$GOOGLE_BUCKET_NAME/work \
  --publish_dir=gs://$GOOGLE_BUCKET_NAME/results \
  --store_dir=gs://$GOOGLE_BUCKET_NAME/store \
  -resume \
  -r main \
  --metadata_tsv https://raw.githubusercontent.com/seandavi/curatedMetagenomicsNextflow/main/samplesheet.test.tsv
```

To view results:

```sh
gsutil ls -larh gs://$GOOGLE_BUCKET_NAME/results
```

To view an individual file:

```sh
gsutil cat PATH_TO_GOOGLE_OBJECT # from above list
```

To cleanup:

```sh
gsutil -m rm -r gs://$GOOGLE_BUCKET_NAME
```

## nf-core tools integration

```sh
pip install nf-core
```

```sh
nf-core launch seandavi/curatedMetagenomicsNextflow
```
