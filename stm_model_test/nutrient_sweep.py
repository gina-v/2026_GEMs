#!/usr/bin/env python
# coding: utf-8

# # Running Nutrient Gradient Sweep Function for Infection Medium Tuning
# 
# ```
# see Obsidian Notes:
# 08.07.26, 08.08.26, 08.10.26, 08.11.26
# 
# continued from:
# /home/gmvaz/2026_GEMs/stm_model_test/080726_081426.ipynb
# 
# tmux: gems
# conda environment: gems
# ```

# # Generate Constrained Models of the STM only condition and the Mixed STM condition


import sys

model_name = sys.argv[1]
medium_name = sys.argv[2]
nutrient_rxn = sys.argv[3]
output_file = sys.argv[4]

# Set the gurobi solver 

import os
from cobra.core.configuration import Configuration

os.environ["GRB_LICENSE_FILE"] = "/home/gmvaz/.gurobi.lic"
Configuration().solver = "gurobi"

# Read in the gene weights CSV files

import pandas as pd

stm_threshold_10_90 = pd.read_csv("../imat_prep_test/stm_threshold_10_90.csv")
stm_threshold_15_85 = pd.read_csv("../imat_prep_test/stm_threshold_15_85.csv")
stm_threshold_25_75 = pd.read_csv("../imat_prep_test/stm_threshold_25_75.csv")

mixed_stm_threshold_10_90 = pd.read_csv("../imat_prep_test/mixed_stm_threshold_10_90.csv")
mixed_stm_threshold_15_85 = pd.read_csv("../imat_prep_test/mixed_stm_threshold_15_85.csv")
mixed_stm_threshold_25_75 = pd.read_csv("../imat_prep_test/mixed_stm_threshold_25_75.csv")

# Make a python set of each gene weights CSV

stm_set_10_90 = set(stm_threshold_10_90["genes"].astype(str).str.strip())
stm_set_15_85 = set(stm_threshold_15_85["genes"].astype(str).str.strip())
stm_set_25_75 = set(stm_threshold_25_75["genes"].astype(str).str.strip())

mixed_stm_set_10_90 = set(mixed_stm_threshold_10_90["genes"].astype(str).str.strip())
mixed_stm_set_15_85 = set(mixed_stm_threshold_15_85["genes"].astype(str).str.strip())
mixed_stm_set_25_75 = set(mixed_stm_threshold_25_75["genes"].astype(str).str.strip())

# QC: Check the count of genes

print("STM1090 gene count",len(stm_set_10_90), "|", "STM1585 gene count", len(stm_set_15_85), "|", "STM2575 gene count", len(stm_set_25_75))
print("Mixed STM1090 gene count",len(mixed_stm_set_10_90), "|", "Mixed STM1585 gene count", len(mixed_stm_set_15_85), "|", "Mixed STM2575 gene count", len(mixed_stm_set_25_75))


# In[2]:


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


# In[3]:


# Generate constrained models with iMAT by integrating the reaction weights

from imatpy.model_creation import generate_model

stm1090_model = generate_model(model = control_model, rxn_weights = rxn_weights_stm1090, method = "imat_restrictions")
stm1585_model = generate_model(model = control_model, rxn_weights = rxn_weights_stm1585, method = "imat_restrictions")
stm2575_model = generate_model(model = control_model, rxn_weights = rxn_weights_stm2575, method = "imat_restrictions")

mixedstm1090_model = generate_model(model = control_model, rxn_weights = rxn_weights_mixedstm1090, method = "imat_restrictions")
mixedstm1585_model = generate_model(model = control_model, rxn_weights = rxn_weights_mixedstm1585, method = "imat_restrictions")
mixedstm2575_model = generate_model(model = control_model, rxn_weights = rxn_weights_mixedstm2575, method = "imat_restrictions")


# In[4]:


stm1090_model


# # Make the `original_models` dicitonary of models

# In[5]:


# Make a dictionary of the models 

original_models = {
    "stm_1090": stm1090_model,
    "stm_1585": stm1585_model,
    "stm_2575": stm2575_model,
    "mixed_stm1090": mixedstm1090_model,
    "mixed_stm1585": mixedstm1585_model,
    "mixed_stm2575": mixedstm2575_model
}

# Make a dictionary with copies of each original model

fba_models = {
    model_name: model.copy()
    for model_name, model in original_models.items()
}


# In[6]:


fba_models


# # Define the Infection Medium Conditions

# In[7]:


# Make a copy of the default medium from any of the models (stm1090 used)

model_copy = fba_models.get("stm_1090")
default_medium = model_copy.medium

# Make 7 identical infection mediums with oxygen ranging from 0 - 3 at 0.5 increments

A_infection_medium = dict(default_medium)
A_infection_medium.update({
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

B_infection_medium = dict(default_medium)
B_infection_medium.update({
    "EX_o2_e": 0.5, # changing in default medium
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

C_infection_medium = dict(default_medium)
C_infection_medium.update({
    "EX_o2_e": 1.0, # changing in the default medium
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

D_infection_medium = dict(default_medium)
D_infection_medium.update({
    "EX_o2_e": 1.5, # changing from the default medium
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

E_infection_medium = dict(default_medium)
E_infection_medium.update({
    "EX_o2_e": 2.0, # changing from the default medium
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

F_infection_medium = dict(default_medium)
F_infection_medium.update({
    "EX_o2_e": 2.5, # changing from the default medium
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

G_infection_medium = dict(default_medium)
G_infection_medium.update({
    "EX_o2_e": 3.0, # changing from the default medium
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


# # Make `base_medias` - dicitonary of the infection medium's

# In[8]:


base_medias = {
	"A_infection_medium": A_infection_medium,
	"B_infection_medium": B_infection_medium,
	"C_infection_medium": C_infection_medium,
	"D_infection_medium": D_infection_medium,
	"E_infection_medium": E_infection_medium,
	"F_infection_medium": F_infection_medium,
	"G_infection_medium": G_infection_medium
}


# # Make the `nutrient_rxn_list` -  list of reaction_ID's who's value will be changed

# In[9]:


nutrient_rxn_list = list(A_infection_medium.keys())
nutrient_rxn_list.remove('EX_o2_e')
nutrient_rxn_list


# # Make the `nutrient_step_count_list`

# In[10]:


import numpy as np

nutrient_step_count = np.logspace(0.5, 7, num = 15, base = 2.0)

nutrient_step_count_array = np.insert(nutrient_step_count, 0, 0) # add values 0.5 to list of nutrient step count

nutrient_step_count_list = nutrient_step_count_array.tolist()

nutrient_step_count_list


# # Define `run_medium_sweep()`

# In[11]:


def run_medium_sweep(*, model_name, medium_name, nutrient_rxn, nutrient_step_count): 

    results = [ ]

    for nutrient_value in nutrient_step_count:

        model = fba_models[model_name].copy()
        medium = base_medias[medium_name].copy()
        medium[nutrient_rxn] = float(nutrient_value)

        try:
            with model:
                model.medium = medium
                solution = model.optimize()

                if solution.status == "optimal": 
                    objective_value = solution.objective_value
                else:
                    objective_value = np.nan # NaN - Not a Number

                results.append({
                    "model": model_name,
                    "medium_condition": medium_name,
                    "nutrient_rxn": nutrient_rxn,
                    "nutrient_rxn_value": float(nutrient_value),
                    "status": solution.status,
                    "objective_value": objective_value
                })

        except Exception as error:
            results.append({
                "model": model_name,
                "medium_condition": medium_name,
                "nutrient_rxn": nutrient_rxn,
                "nutrient":float(nutrient_value),
                "status": "failed",
                "objective_value": np.nan,
                "error": str(error)
            })

    results_df = pd.DataFrame(results)

    return results_df


# In[12]:


stm1090_ca2_E_medium_results = run_medium_sweep(model_name = "stm_1090", medium_name = "E_infection_medium", nutrient_rxn = "EX_ca2_e", nutrient_step_count = nutrient_step_count_list)
stm1090_ca2_E_medium_results.to_csv("nutrient_sweep_081226/stm1090_ca2_E_medium_results.csv")


# # Define `run_medium_sweep_fast()` with `model.slim_optimize()` instead of `model.optimize`

# In[13]:


def run_medium_sweep_fast(*, model_name, medium_name, nutrient_rxn, nutrient_step_count): 

    results = [ ]

    for nutrient_value in nutrient_step_count:

        model = fba_models[model_name].copy()
        medium = base_medias[medium_name].copy()
        medium[nutrient_rxn] = float(nutrient_value)

        try:
            with model:
                model.medium = medium
                objective_value = model.slim_optimize()

                if isinstance(objective_value, float) == True:
                        status = "optimal"
                else:
                    status = "infeasible"

                results.append({
                    "model": model_name,
                    "medium_condition": medium_name,
                    "nutrient_rxn": nutrient_rxn,
                    "nutrient_rxn_value": float(nutrient_value),
                    "status": status,
                    "objective_value": objective_value
                })

        except Exception as error:
            #raise
            results.append({
                "model": model_name,
                "medium_condition": medium_name,
                "nutrient_rxn": nutrient_rxn,
                "nutrient":float(nutrient_value),
                "status": "failed",
                "objective_value": np.nan,
                "error": str(error)
            })

    results_df = pd.DataFrame(results)

    return results_df


# In[14]:


#stm1090_cbl1_A_medium_results = run_medium_sweep_fast(model_name = "stm_1090", medium_name = "A_infection_medium", nutrient_rxn = "EX_cbl1_e", nutrient_step_count = nutrient_step_count_list)
#stm1090_cbl1_A_medium_results.to_csv("nutrient_sweep_081226/stm1090_cbl1_A_medium_results.csv")

# stm_1090 A_infection_medium EX_cbl1_e nutrient_sweep_081226/stm1090_cbl1_test_medium_results.csv

medium_results = run_medium_sweep_fast(model_name = model_name, medium_name = medium_name, nutrient_rxn = nutrient_rxn, nutrient_step_count = nutrient_step_count_list)

medium_results.to_csv(output_file)

# stm_1090 B_infection_medium EX_cbl1_e nutrient_sweep_081226/stm1090_cbl1_B_medium_results.csv

# stm_1090 C_infection_medium EX_cbl1_e nutrient_sweep_081226/stm1090_cbl1_C_medium_results.csv

# stm_1090 D_infection_medium EX_cbl1_e nutrient_sweep_081226/stm1090_cbl1_D_medium_results.csv

# stm_1090 E_infection_medium EX_cbl1_e nutrient_sweep_081226/stm1090_cbl1_E_medium_results.csv























