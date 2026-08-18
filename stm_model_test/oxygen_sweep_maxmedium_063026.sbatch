#!/bin/bash
#SBATCH --job-name=oxygen_sweep_063026
#SBATCH --partition=bml
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/home/gmvaz/2026_GEMs/stm_model_test/oxygen_sweep_results/063026/oxygen_sweep_%j.out
#SBATCH --error=/home/gmvaz/2026_GEMs/stm_model_test/oxygen_sweep_results/063026/oxygen_sweep_%j.err

set -euo pipefail

echo "Job started on: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"

# Go to project directory
cd /home/gmvaz/2026_GEMs/stm_model_test

# Load your normal shell setup, including conda setup
source ~/.bashrc

# Activate conda environment
conda activate gems

# Tell Gurobi where the license file is
export GRB_LICENSE_FILE="$HOME/.gurobi.lic"

# Variables used by the Python script
export STM_MODEL_TEST_DIR=/home/gmvaz/2026_GEMs/stm_model_test
export COBRA_SOLVER=gurobi
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Run script
python -u oxygen_sweep_maxmedium_063026.py

echo "End time: $(date)"
echo "Job finished."
