#!/bin/bash
#SBATCH --job-name=context_models_031526
#SBATCH --output=/home/gmvaz/2026_GEMs/stm_model_test/run_031526/logs/context_models_%j.out
#SBATCH --error=/home/gmvaz/2026_GEMs/stm_model_test/run_031526/logs/context_models_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=bmh

source ~/.bashrc
conda activate gems

export GRB_LICENSE_FILE=$HOME/.gurobi.lic

cd /home/gmvaz/2026_GEMs/stm_model_test

python generate_models_0315.py

#python gen_models2.py --threads 16 /home/gmvaz/2026_GEMs/imat_prep_test/stm_threshold_10_90.csv
