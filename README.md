# README

google.api_core.exceptions.PermissionDenied: 403 Permission monitoring.metricDescriptors.create denied (or the resource may not exist).

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
  - `NCBI_accession`, a semicolon-separated list of SRRs
- Can be a file or a web url

If using a Google Bucket, the name bucket must not have underscores.

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
