#!/bin/bash -l
# submit_alpine.sh
#
# Usage:
#   sbatch submit_alpine.sh <RUN_IDS> <SAMPLE_ID> [REVISION]
#
# Examples:
#   sbatch submit_alpine.sh SRR1234567 sample_001
#   sbatch submit_alpine.sh SRR1234567 sample_001 refactor/pipeline-clarity-resilience
#
# The optional REVISION is passed to `nextflow run -r` and may be a branch,
# tag, or commit SHA. This is the preferred way to test branch-specific changes
# on the cluster without staging a local working copy into scratch space.

#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal

set -euo pipefail
set -x

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: sbatch submit_alpine.sh <RUN_IDS> <SAMPLE_ID> [REVISION]" >&2
    exit 1
fi

RUN_IDS="$1"
SAMPLE_ID="$2"
REVISION="${3:-main}"

export GOOGLE_APPLICATION_CREDENTIALS="$HOME/curatedmetagenomicdata-232f4a306d1d.json"

module load singularity
module load git
module load nextflow

if [[ -n "${SLURM_SCRATCH:-}" ]]; then
    cd "$SLURM_SCRATCH"
fi

export NXF_MODE=google

nextflow run seandavi/curatedMetagenomicsNextflow \
    -r "$REVISION" \
    --run_ids "$RUN_IDS" \
    --sample_id "$SAMPLE_ID" \
    -profile alpine
