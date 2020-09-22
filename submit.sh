#!/bin/bash
module load singularity
module load git
export GOOGLE_APPLICATION_CREDENTIALS=/data/sdavis2/projects/curatedMetagenomicsNextflow/curatedmetagenomicdata-f056e9ef05f6.json
WORKDIR=nf_`uuidgen`

echo "using WORKDIR $WORKDIR"

mkdir $WORKDIR
cd $WORKDIR
cp ../main.nf .
cp ../nextflow.config .
export NXF_MODE=google
RUNS=$1 # semicolon-separated
SAMP=$2 # just name
STUDY=$3
NXF_MODE=google nextflow run main.nf -with-weblog https://cmgd-telemetry-whnnxetv4q-uc.a.run.app/nextflow/events -profile nih_slurm --runs $RUNS --samp $SAMP --study $STUDY
cd ..
rm -rf $WORKDIR
