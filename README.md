# README

## The pipeline, in brief

- fastq-dump
- setup databases
  - metaphlan
  - chocophlan
  - uniref
- metaphlan bug list
- metaphlan markers
- humann (including various aggregations)

## Inputs and Outputs

The `metadata_tsv` file must be:

- tab-separated
- must contain columns
  - `sample_id`
  - `study_name`
  - `NCBI_accessions`, a semicolon-separated list of SRRs
- Can be a file or a web url

## On sample ids

We use a simple approach to create sample ids. The `study_name` and `sample_id` are
first concatenated by `::`. Then, we base64 encode. For example:

```sh
echo 'study_name1::sample_name1' | base64
```

This yields:

```
c3R1ZHlfbmFtZTE6OnNhbXBsZV9uYW1lMQo=
```

To decode a sample id:

```sh
echo 'c3R1ZHlfbmFtZTE6OnNhbXBsZV9uYW1lMQo=' | base64 -d
```

which gives back the original string:

```
study_name1::sample_name1
```




## Install

```sh
export NXF_MODE=google
curl https://get.nextflow.io | bash
```

## Google setup

You will need to be able to access google cloud storage as well as to 
run the Google Cloud Pipeline API. This requires credentials to do so.
You can either use Google Default Application Credentials or a key file.
The latter is the recommended approach. If you need a keyfile, contact
the person who owns the Google Project you'll be using. 

### Keyfile setup

Once you get a keyfile (which is a json file), set the environment variable:

```sh
export GOOGLE_APPLICATION_CREDENTIALS=/path/to/keyfile/keyfile-name.json
```

## Execution

This assumes that you are running on Google, that credentials are set up correctly,
and that you have a Google Storage Bucket already created.

```sh
export GOOGLE_BUCKET_NAME="YOUR_GOOGLE_BUCKET"
# optional
gsutil mb gs://$GOOGLE_BUCKET_NAME
```

You can now run test data. 

```sh
nextflow run seandavi/curatedMetagenomicsNextflow \
  -profile google \
  -work-dir gs://$GOOGLE_BUCKET_NAME/work \
  --publish_dir=gs://$GOOGLE_BUCKET_NAME/results \
  --store_dir=gs://$GOOGLE_BUCKET_NAME/store \
  -resume \
  --metadata_tsv https://raw.githubusercontent.com/seandavi/curatedMetagenomicsNextflow/main/samplesheet.test.tsv
```

## Viewing results

```sh
gsutil ls -lahr gs://$GOOGLE_BUCKET_NAME/results
```
