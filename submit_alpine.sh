#!/bin/bash
#SBATCH --mem=48G
#SBATCH --cpus-per-task=16
set -x

echo "working in $SLURM_SCRATCH"

cp main.nf $SLURM_SCRATCH
cp nextflow.config $SLURM_SCRATCH


export GOOGLE_APPLICATION_CREDENTIALS=$HOME/omicidx-338300-cbd1527c319e.json
module load singularity
module load git
module load nextflow

cd $SLURM_SCRATCH
export NXF_MODE=google
nextflow run main.nf --run_ids=$1 --sample_id=$2 -profile alpine
