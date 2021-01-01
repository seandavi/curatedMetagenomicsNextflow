#!/bin/bash
set -x 

module load singularity
module load git
export GOOGLE_APPLICATION_CREDENTIALS=/data/sdavis2/projects/curatedMetagenomicsNextflow/curatedmetagenomicdata-f056e9ef05f6.json
WORKDIR=/lscratch/$SLURM_JOB_ID/

echo "using WORKDIR $WORKDIR"

cp main.nf $WORKDIR
cp nextflow.config $WORKDIR
cp nextflow $WORKDIR
cd $WORKDIR
export NXF_MODE=google
#NXF_MODE=google nextflow run main.nf -with-weblog https://cmgd-telemetry-whnnxetv4q-uc.a.run.app/nextflow/events -profile nih_slurm --metadata_tsv $1 && cd .. && rm -rf $WORKDIR
NXF_MODE=google nextflow run main.nf -with-weblog https://cmgd-telemetry-whnnxetv4q-uc.a.run.app/nextflow/events -profile nih_node_only --metadata_tsv --srr $1 --uuid $2 
