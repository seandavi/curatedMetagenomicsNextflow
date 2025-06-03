#!/bin/bash
# submit_unitn.sh

# Usage
# qsub -N job_name -v metadata_tsv=/absolute/path/to/samples.tsv submit_unitn.sh

# 1 node
# 1 core
# 24 GB
#PBS -l select=1:ncpus=1:mem=24gb

# execution queue: common_cpuQ; CIBIO_cpuQ
#PBS -q CIBIO_cpuQ

# name the job on the command line: $ qsub -N job_name submit_unitn.sh

# write output and error as job is progressing
#PBS -k oed

# send email on job abort, begin, end
##PBS -m abe

# email address for notifications
##PBS -M <your_email_address>

set -x

echo "working in $UNITN_SCRATCH"

echo "metadata_tsv: $metadata_tsv"

export GOOGLE_APPLICATION_CREDENTIALS=/absolute/path/to/keyfile.json

# for allowing singularity to access $HOME/.ncbi/user-settings.mkfg
# still not sure why user-settings that is supposed to be in Docker container isn't accessible
export NXF_SINGULARITY_HOME_MOUNT=true

# activate conda environment
conda activate metagenomicsMAC

# load singularity
module load singularity-3.4.0

cd $UNITN_SCRATCH
export NXF_MODE=google

nextflow run ASAP-MAC/metagenomicsNextflowMAC --metadata_tsv=$metadata_tsv -profile unitn -with-weblog https://nf-telemetry-819875667022.us-central1.run.app/nextflow-telemetry/events
