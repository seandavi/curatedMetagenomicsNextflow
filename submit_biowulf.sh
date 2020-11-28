#!/bin/bash
module load singularity
module load git
export GOOGLE_APPLICATION_CREDENTIALS=/data/sdavis2/projects/curatedMetagenomicsNextflow/curatedmetagenomicdata-f056e9ef05f6.json
WORKDIR=nf_`uuidgen`

echo "using WORKDIR $WORKDIR"

mkdir $WORKDIR
cd $WORKDIR
cp ../$1 .
cp ../main.nf .
cp ../nextflow.config .
export NXF_MODE=google
NXF_MODE=google nextflow run main.nf -with-weblog https://cmgd-telemetry-whnnxetv4q-uc.a.run.app/nextflow/events -profile nih_slurm --metadata_tsv $1 && cd .. && rm -rf $WORKDIR
