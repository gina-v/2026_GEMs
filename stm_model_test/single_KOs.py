import sys
model_name = sys.argv[1]

import os
from cobra.core.configuration import Configuration

# path to gurobi license file
os.environ["GRB_LICENSE_FILE"] = "/home/gmvaz/.gurobi.lic"

# change default cobra solver to gurobi
Configuration().solver = "gurobi"

# Read in the gene weights CSV files
import pandas as pd
stm_threshold_10_90 = pd.read_csv("../imat_prep_test/stm_threshold_10_90.csv")
stm_threshold_15_85 = pd.read_csv("../imat_prep_test/stm_threshold_15_85.csv")
stm_threshold_25_75 = pd.read_csv("../imat_prep_test/stm_threshold_25_75.csv")

mixed_stm_threshold_10_90 = pd.read_csv("../imat_prep_test/mixed_stm_threshold_10_90.csv")
mixed_stm_threshold_15_85 = pd.read_csv("../imat_prep_test/mixed_stm_threshold_15_85.csv")
mixed_stm_threshold_25_75 = pd.read_csv("../imat_prep_test/mixed_stm_threshold_25_75.csv")

# Make a set of each gene weights CSV
stm_set_10_90 = set(stm_threshold_10_90["genes"].astype(str).str.strip())
stm_set_15_85 = set(stm_threshold_15_85["genes"].astype(str).str.strip())
stm_set_25_75 = set(stm_threshold_25_75["genes"].astype(str).str.strip())

mixed_stm_set_10_90 = set(mixed_stm_threshold_10_90["genes"].astype(str).str.strip())
mixed_stm_set_15_85 = set(mixed_stm_threshold_15_85["genes"].astype(str).str.strip())
mixed_stm_set_25_75 = set(mixed_stm_threshold_25_75["genes"].astype(str).str.strip())

# QC: Check the count of genes
print("STM1090 gene count",len(stm_set_10_90), "|", "STM1585 gene count", len(stm_set_15_85), "|", "STM2575 gene count", len(stm_set_25_75))
print("Mixed STM1090 gene count",len(mixed_stm_set_10_90), "|", "Mixed STM1585 gene count", len(mixed_stm_set_15_85), "|", "Mixed STM2575 gene count", len(mixed_stm_set_25_75))

# Load STM_V1_0
from cobra.io import load_model, load_json_model

control_model = load_model("STM_v1_0")

# Save all genes in the control model in a set
model_genes = set([g.id for g in control_model.genes])

# QC: Check the count of genes
print("STM_V1_0 Control Model has", len(model_genes), "total genes")

# Subset the gene weights lists for each sample to keep only the genes also found in the control model
stm1090_model = stm_threshold_10_90[stm_threshold_10_90["genes"].isin(model_genes)]
stm1585_model = stm_threshold_15_85[stm_threshold_15_85["genes"].isin(model_genes)]
stm2575_model = stm_threshold_25_75[stm_threshold_25_75["genes"].isin(model_genes)]

mixed_stm1090_model = mixed_stm_threshold_10_90[mixed_stm_threshold_10_90["genes"].isin(model_genes)]
mixed_stm1585_model = mixed_stm_threshold_15_85[mixed_stm_threshold_15_85["genes"].isin(model_genes)]
mixed_stm2575_model = mixed_stm_threshold_25_75[mixed_stm_threshold_25_75["genes"].isin(model_genes)]

# Align the gene weights lists to the control model genes list and assign missing genes a weight of 0
subset_weights_stm1090 = stm1090_model.set_index("genes")["threshold_category"].reindex(model_genes, fill_value=0)
subset_weights_stm1585 = stm1585_model.set_index("genes")["threshold_category"].reindex(model_genes, fill_value=0)
subset_weights_stm2575 = stm2575_model.set_index("genes")["threshold_category"].reindex(model_genes, fill_value=0)

subset_weights_mixedstm1090 = mixed_stm1090_model.set_index("genes")["threshold_category"].reindex(model_genes, fill_value=0)
subset_weights_mixedstm1585 = mixed_stm1585_model.set_index("genes")["threshold_category"].reindex(model_genes, fill_value=0)
subset_weights_mixedstm2575 = mixed_stm2575_model.set_index("genes")["threshold_category"].reindex(model_genes, fill_value=0)


# QC: Confirm all samples have the same amount of gene weights
# 1. Make a list of the subset gene weights
subset_gene_weights = [
    subset_weights_stm1090,
    subset_weights_stm1585,
    subset_weights_stm2575,
    subset_weights_mixedstm1090,
    subset_weights_mixedstm1090,
    subset_weights_mixedstm2575]

# 2. Make a list of how many genes are in each list
gene_counts = [len(weights) for weights in subset_gene_weights]

# 3. Compare how many genes are in each list and print the total if they are all the same
if len(set(gene_counts)) == 1:
    print(f"All models have {gene_counts[0]} genes")
else:
    print("Gene counts differ between models")
    print(gene_counts)

# Convert gene weights to reaction weights with iMAT

from imatpy.parse_gpr import gene_to_rxn_weights
rxn_weights_stm1090 = gene_to_rxn_weights(control_model, subset_weights_stm1090)
rxn_weights_stm1585 = gene_to_rxn_weights(control_model, subset_weights_stm1585)
rxn_weights_stm2575 = gene_to_rxn_weights(control_model, subset_weights_stm2575)

rxn_weights_mixedstm1090 = gene_to_rxn_weights(control_model, subset_weights_mixedstm1090)
rxn_weights_mixedstm1585 = gene_to_rxn_weights(control_model, subset_weights_mixedstm1585)
rxn_weights_mixedstm2575 = gene_to_rxn_weights(control_model, subset_weights_mixedstm2575)

# QC: Confirm all samples have the same amount of reaction weights
# 1. Make a list of the reaction weights
reaction_weights = [
    rxn_weights_stm1090,
    rxn_weights_stm1585,
    rxn_weights_stm2575,
    rxn_weights_mixedstm1090,
    rxn_weights_mixedstm1585,
    rxn_weights_mixedstm2575]

# 2. Make a list of how many genes are in each list
reactions_counts = [len(weights) for weights in reaction_weights]

# 3. Compare how many genes are in each list and print the total if they are all the same
if len(set(reactions_counts)) == 1:
    print(f"All models have {reactions_counts[0]} reactions")
else:
    print("Reaction counts differ between models")
    print(reactions_counts)

# Generate constrained models with iMAT by integrating the reaction weights
from imatpy.model_creation import generate_model

rxn_weights = None
if model_name=="stm_1090":
    rxn_weights = rxn_weights_stm1090
elif model_name == "stm1585":
    rxn_weights = rxn_weights_stm1585
elif model_name == "stm2575":
    rxn_weights = rxn_weights_stm2575
elif model_name == "mixed_stm1090":
    rxn_weights = rxn_weights_mixedstm1090
elif model_name == "mixed_stm1585":
    rxn_weights = rxn_weights_mixedstm1585
elif model_name == "mixed_stm2575":
    rxn_weights = rxn_weights_mixedstm2575
elif model_name == "stm_unconstrained":
    rxn_weights = None

assert rxn_weights is not None

model = generate_model(model = control_model, rxn_weights = rxn_weights, method = "imat_restrictions")

# make a dictionary of all the models generated 
models = {
    model_name: model
}

# Make a copy of the default medium from stm_v1_0
default_medium = dict(models["stm_unconstrained"].medium)

# Add mouse-derived and host-derived gut metabolites to the default medium to create a mouse gut infection medium
infection_medium = dict(default_medium)
infection_medium.update({
    "EX_o2_e": 0.0, # removing from default medium
    "EX_tet_e": 4.0,
    "EX_no3_e": 0.5,
    "EX_etha_e": 1.0,
    "EX_12ppd__R_e": 0.1,
    "EX_12ppd__S_e": 0.1,
    "EX_lac__L_e": 10.0,
    "EX_for_e": 4.0,
    "EX_succ_e": 3.0,
    "EX_asp__L_e": 5.0,
    "EX_mal__L_e": 0.2,
    "EX_galct__D_e": 10.0,
    "EX_glcr_e": 10.0,
    "EX_pyr_e": 10.0,
    "EX_arab__L_e": 13.3,
    "EX_but_e": 10.0,
    "EX_ac_e": 50.0,
    "EX_ppa_e": 5.0,
    "EX_glc__D_e": 0.0 # removing from default medium
})

import pandas as pd 
from time import time 

from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion)

geneko_results_list = [ ]

for name, model in models.items():

    geneko_model = model.copy()
    geneko_medium = infection_medium.copy()
    geneko_model.medium = geneko_medium

    geneko_results = single_gene_deletion(geneko_model, geneko_model.genes, method="fba")

    geneko_results['model'] = name
    geneko_results_list.append(geneko_results)
    geneko_results_df = pd.concat(geneko_results_list)
    geneko_results_df.to_csv(f"geneko_results.csv", index=False) 


rxnko_results_list = [ ]

for name, model in models.items():

    rxnko_model = model.copy()
    rxnko_medium = infection_medium.copy()
    rxnko_model.medium = rxnko_medium

    rxnko_results = single_reaction_deletion(rxnko_model, rxnko_model.reactions, method="fba")

    rxnko_results['model'] = name
    rxnko_results_list.append(rxnko_results)
    rxnko_results_df = pd.concat(rxnko_results_list)
    rxnko_results_df.to_csv(f"rxnko_results.csv", index=False)

# add the rxn name and gene name to the dataframe columns
# check that the for loops work
# check if I need system args
# if I need them then how do they work?? 
# how to run with sbatch and snakemake
                  