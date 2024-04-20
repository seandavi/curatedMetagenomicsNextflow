#!/bin/bash
#SBATCH --mem=48G
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
set -x

NXF_SINGULARITY_CACHEDIR=$SCRATCH/singularity_cache
mkdir -p $NXF_SINGULARITY_CACHEDIR

PARENT_DIR=$SCRATCH/cmgd_data
WORKDIR=$PARENT_DIR/$SLURM_JOB_ID
mkdir -p $WORKDIR

echo "working in $WORKDIR"

cp main.nf $WORKDIR
cp nextflow.config $WORKDIR


export GOOGLE_APPLICATION_CREDENTIALS=$HOME/omicidx-338300-cbd1527c319e.json
module load nextflow

cd $WORKDIR
export NXF_MODE=google
nextflow run main.nf --run_ids=$1 --sample_id=$2 -profile anvil

cd $HOME

# cleanup
#rm -r $WORKDIR
