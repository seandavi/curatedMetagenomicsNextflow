#!/bin/bash
set -x 

module load singularity
module load git
export GOOGLE_APPLICATION_CREDENTIALS=/home/seandavi/curatedmetagenomicdata-f056e9ef05f6.jsonexport 

WORKDIR=$LOCAL/

echo "using WORKDIR $WORKDIR"

cp main.nf $WORKDIR
cp nextflow.config $WORKDIR
cp nextflow $WORKDIR
cd $WORKDIR
export NXF_MODE=google
#NXF_MODE=google nextflow run main.nf -with-weblog https://cmgd-telemetry-whnnxetv4q-uc.a.run.app/nextflow/events -profile nih_slurm --metadata_tsv $1 && cd .. && rm -rf $WORKDIR
NXF_MODE=google ./nextflow run main.nf -with-weblog https://cmgd-telemetry-whnnxetv4q-uc.a.run.app/nextflow/events -profile bridges_node_only --srr $1 --uuid $2 
