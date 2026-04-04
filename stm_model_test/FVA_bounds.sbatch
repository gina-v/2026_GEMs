#!/bin/bash

#SBATCH --job-name=context_models_031726
#SBATCH --output=/home/gmvaz/2026_GEMs/stm_model_test/run_031526/logs/context_models_%j.out
#SBATCH --error=/home/gmvaz/2026_GEMs/stm_model_test/run_031526/logs/context_models_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=bmh

source ~/.bashrc
conda activate gems

export GRB_LICENSE_FILE=$HOME/.gurobi.lic

cd /home/gmvaz/2026_GEMs/stm_model_test

#03/23/26 - FVA
python generate_models3.py /home/gmvaz/2026_GEMs/imat_prep_test/stm_threshold_10_90.csv FVA --output run_031526/context_models/stm1090_FVA.json
#python generate_models3.py /home/gmvaz/2026_GEMs/imat_prep_test/stm_threshold_15_85.csv FVA --output run_031526/context_models/stm1585_FVA.json
#python generate_models3.py /home/gmvaz/2026_GEMs/imat_prep_test/stm_threshold_25_75.csv FVA --output run_031526/context_models/stm2575_FVA.json
#python generate_models3.py /home/gmvaz/2026_GEMs/imat_prep_test/mixed_stm_threshold_10_90.csv FVA --output run_031526/context_models/mixed1090_FVA.json
#python generate_models3.py /home/gmvaz/2026_GEMs/imat_prep_test/mixed_stm_threshold_15_85.csv FVA --output run_031526/context_models/mixed1585_FVA.json
#python generate_models3.py /home/gmvaz/2026_GEMs/imat_prep_test/mixed_stm_threshold_25_75.csv FVA --output run_031526/context_models/mixed2575_FVA.json

#03/23/26 - MILP
#python generate_models3.py /home/gmvaz/2026_GEMs/imat_prep_test/stm_threshold_10_90.csv MILP --output run_031526/context_models/stm1090_MILP.json
#python generate_models3.py /home/gmvaz/2026_GEMs/imat_prep_test/stm_threshold_15_85.csv MILP --output run_031526/context_models/stm1585_MILP.json
#python generate_models3.py /home/gmvaz/2026_GEMs/imat_prep_test/stm_threshold_25_75.csv MILP --output run_031526/context_models/stm2575_MILP.json
#python generate_models3.py /home/gmvaz/2026_GEMs/imat_prep_test/mixed_stm_threshold_10_90.csv MILP --output run_031526/context_models/mixed1090_MILP.json
#python generate_models3.py /home/gmvaz/2026_GEMs/imat_prep_test/mixed_stm_threshold_15_85.csv MILP --output run_031526/context_models/mixed1585_MILP.json
#python generate_models3.py /home/gmvaz/2026_GEMs/imat_prep_test/mixed_stm_threshold_25_75.csv MILP --output run_031526/context_models/mixed2575_MILP.json


