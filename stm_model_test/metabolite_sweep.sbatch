#!/bin/bash
#SBATCH --job-name=metabolite_sweep_test
#SBATCH --output=/home/gmvaz/2026_GEMs/stm_model_test/logs/metabolite_sweep_%j.out
#SBATCH --error=/home/gmvaz/2026_GEMs/stm_model_test/logs/metabolite_sweep_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --partition=bmh

source ~/.bashrc
conda activate gems

snakemake -s Snakefile.py -j 16  