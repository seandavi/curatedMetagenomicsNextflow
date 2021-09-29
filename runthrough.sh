
# install nextflow
export NXF_MODE=google
curl https://get.nextflow.io | bash

# gsutil if you like
# pip install gsutil

# setup google access -- YOU NEED TO EDIT
export GOOGLE_APPLICATION_CREDENTIALS=/Users/seandavis/Documents/git/curatedMetagenomicsNextflow/curatedmetagenomicdata-0ca253fe73b7.json

# where to put stuff -- YOU NEED TO EDIT
export GOOGLE_BUCKET_NAME=temp-testing-cmgd
# make bucket if it doesn't exist
gsutil mb gs://$GOOGLE_BUCKET_NAME

# pull the workflow directly from github
nextflow pull seandavi/curatedMetagenomicsNextflow -r main # could be a specific revision

# And run with sample data
nextflow run seandavi/curatedMetagenomicsNextflow \
  -profile google \
  -work-dir gs://$GOOGLE_BUCKET_NAME/work \
  -r main \
  --metadata_tsv https://raw.githubusercontent.com/seandavi/curatedMetagenomicsNextflow/main/samplesheet.test.tsv \
  --bucket gs://$GOOGLE_BUCKET_NAME