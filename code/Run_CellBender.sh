#!/bin/bash -l

#SBATCH --job-name=cellbender
#SBATCH --output=logs/cellbender_%A_%a.out
#SBATCH --error=logs/cellbender_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=46G
#SBATCH --array=0-23
#SBATCH --gres=gpu:1
#SBATCH --partition=l40s

# Run from the repo root, or set REPO_DIR explicitly.
REPO_DIR="${REPO_DIR:-$(git rev-parse --show-toplevel)}"
RAW_DIR="${REPO_DIR}/data/raw_h5"
OUT_DIR="${REPO_DIR}/data/cellbender"

conda activate cellbender

names=($(cd "$RAW_DIR" && ls *.h5 | sed 's/\.h5$//'))
name=${names[${SLURM_ARRAY_TASK_ID:-0}]}

mkdir -p "${OUT_DIR}/${name}_cellbender"

cellbender remove-background \
  --input  "${RAW_DIR}/${name}.h5" \
  --output "${OUT_DIR}/${name}_cellbender/${name}_cellbender.h5" \
  --cuda \
  --expected-cells 5000
