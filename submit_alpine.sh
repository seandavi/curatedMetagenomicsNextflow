#!/bin/bash
#SBATCH --mem=48G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00

# Usage:
# sbatch submit_anvil.sh METADATA.TSV_FILE
set -x

echo "working in $SLURM_SCRATCH"

cp main.nf $SLURM_SCRATCH
cp nextflow.config $SLURM_SCRATCH


export GOOGLE_APPLICATION_CREDENTIALS=$HOME/curatedmetagenomicdata-232f4a306d1d.json
module load singularity
module load git
module load nextflow

cd $SLURM_SCRATCH
export NXF_MODE=google
#nextflow run main.nf --run_ids=$1 --sample_id=$2 -profile alpine
nextflow run seandavi/curatedMetagenomicsNextflow --metadata_tsv=$1 -profile alpine -with-weblog https://nf-telemetry-819875667022.us-central1.run.app/nextflow-telemetry/events
