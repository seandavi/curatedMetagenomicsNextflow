#!/bin/bash
# submit_alpine.sh
#
# Usage:
#   sbatch submit_alpine.sh <RUN_IDS> <SAMPLE_ID>
#
# Example:
#   sbatch submit_alpine.sh SRR1234567 sample_001
#
# This script stages the current local pipeline checkout into node-local scratch
# and runs that staged copy. That keeps cluster submissions aligned with the
# branch and local edits being tested, including modularized `modules/` and
# `conf/` directories that are now required by the pipeline entrypoint.

#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00

set -euo pipefail
set -x

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch submit_alpine.sh <RUN_IDS> <SAMPLE_ID>" >&2
    exit 1
fi

RUN_IDS="$1"
SAMPLE_ID="$2"

if [[ -z "${SLURM_SCRATCH:-}" ]]; then
    echo "SLURM_SCRATCH is not set" >&2
    exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="$SLURM_SCRATCH/curatedMetagenomicsNextflow_${SLURM_JOB_ID}"

mkdir -p "$WORKDIR"

echo "staging pipeline from $REPO_ROOT to $WORKDIR"

cp "$REPO_ROOT/main.nf" "$WORKDIR/"
cp "$REPO_ROOT/nextflow.config" "$WORKDIR/"
cp "$REPO_ROOT/submit_alpine.sh" "$WORKDIR/"
cp -R "$REPO_ROOT/conf" "$WORKDIR/"
cp -R "$REPO_ROOT/modules" "$WORKDIR/"

export GOOGLE_APPLICATION_CREDENTIALS="$HOME/curatedmetagenomicdata-232f4a306d1d.json"

module load singularity
module load git
module load nextflow

cd "$WORKDIR"

export NXF_MODE=google

nextflow run main.nf \
    --run_ids "$RUN_IDS" \
    --sample_id "$SAMPLE_ID" \
    -profile alpine
