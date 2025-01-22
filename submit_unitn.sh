#!/bin/bash
# submit_unitn.sh

# Usage
# For a single sample: use line 60
# qsub -N job_name -v run_ids=ERR4330026,sample_id=SAMEA7041133 submit_unitn.sh
# For multiple samples: use line 61
# qsub -N job_name -v metadata_tsv=samples.tsv submit_unitn.sh

# 2 node
# 4 core
# 24 GB
#PBS -l select=2:ncpus=4:mem=24gb

# 8 hours maximum execution time
#PBS -l walltime=08:00:00

# execution queue: common_cpuQ; CIBIO_cpuQ
#PBS -q common_cpuQ

# name the job on the command line: $ qsub -N job_name submit_unitn.sh

# merge stdout and stderr
##PBS -j oe

# write output and error as job is progressing
#PBS -k oed

# send email on job abort, begin, end
#PBS -m abe

# email address for notifications
##PBS -M <your_email_address>

set -x

echo "working in $UNITN_SCRATCH"

#echo "run_ids: $run_ids"
#echo "sample_id: $sample_id"
echo "metadata_tsv: $metadata_tsv"

#cp main.nf $UNITN_SCRATCH
#cp nextflow.config $UNITN_SCRATCH

export GOOGLE_APPLICATION_CREDENTIALS=/home/kaelyn.long/google_cred/curatedmetagenomicdata-232f4a306d1d.json

# for allowing singularity to access $HOME/.ncbi/user-settings.mkfg
# still not sure why user-settings that is supposed to be in Docker container isn't accessible
export NXF_SINGULARITY_HOME_MOUNT=true

# activate conda environment
conda activate metagenomicsMAC

# load singularity
module load singularity-3.4.0

cd $UNITN_SCRATCH
export NXF_MODE=google
#nextflow run main.nf --run_ids=$1 --sample_id=$2 -profile unitn
#nextflow run seandavi/curatedMetagenomicsNextflow --run_ids=$run_ids --sample_id=$sample_id -profile unitn -with-weblog https://nf-telemetry-819875667022.us-central1.run.app/nextflow-telemetry/events
nextflow run seandavi/curatedMetagenomicsNextflow --metadata_tsv=$metadata_tsv -profile unitn -with-weblog https://nf-telemetry-819875667022.us-central1.run.app/nextflow-telemetry/events
