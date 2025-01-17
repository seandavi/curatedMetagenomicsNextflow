#!/bin/bash
#SBATCH --mem=48G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --account=bio240036
set -x

NXF_SINGULARITY_CACHEDIR=$SCRATCH/singularity_cache
mkdir -p $NXF_SINGULARITY_CACHEDIR

PARENT_DIR=$SCRATCH/cmgd_data
WORKDIR=$PARENT_DIR/$SLURM_JOB_ID
mkdir -p $WORKDIR

echo "working in $WORKDIR"

export GOOGLE_APPLICATION_CREDENTIALS=$HOME/curatedmetagenomicdata-232f4a306d1d.json
module load nextflow

cd $WORKDIR
export NXF_MODE=google
nextflow run seandavi/curatedMetagenomicsNextflow --run_ids=$1 --sample_id=$2 -profile anvil -with-weblog https://efc7-71-244-146-119.ngrok-free.app/nextflow-telemetry/events
nextflow clean
