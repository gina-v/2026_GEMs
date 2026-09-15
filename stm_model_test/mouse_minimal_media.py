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

#stm1090_model = generate_model(model = control_model, rxn_weights = rxn_weights_stm1090, method = "imat_restrictions")
#stm1585_model = generate_model(model = control_model, rxn_weights = rxn_weights_stm1585, method = "imat_restrictions")
#stm2575_model = generate_model(model = control_model, rxn_weights = rxn_weights_stm2575, method = "imat_restrictions")

mixedstm1090_model = generate_model(model = control_model, rxn_weights = rxn_weights_mixedstm1090, method = "imat_restrictions")
mixedstm1585_model = generate_model(model = control_model, rxn_weights = rxn_weights_mixedstm1585, method = "imat_restrictions")
mixedstm2575_model = generate_model(model = control_model, rxn_weights = rxn_weights_mixedstm2575, method = "imat_restrictions")

# %%
# make a dictionary of all the models 

models = {
    "stm_unconstrained": control_model,
    "mixed_stm1090": mixedstm1090_model,
    "mixed_stm1585": mixedstm1585_model,
    "mixed_stm2575": mixedstm2575_model,
}

# Make a copy of the default medium from stm_v1_0
default_medium = dict(models["stm_unconstrained"].medium)

# Add mouse-derived, gut metabolites to the default medium to create a nouse gut minimal medium

mouse_min_medium = dict(default_medium)
mouse_min_medium.update({
    "EX_o2_e": 0.5, # updating to anoxic conditions
    "EX_no3_e": 0.5,
    "EX_for_e": 4.0,
    "EX_asp__L_e": 5.0,
    "EX_mal__L_e": 0.2,
    "EX_galct__D_e": 10.0,
    "EX_glcr_e": 10.0,
    "EX_pyr_e": 10.0,
    "EX_arab__L_e": 13.3,
    "EX_glc__D_e": 0.0 # removing from default medium
})

import pandas as pd

results_dir = ("/home/gmvaz/2026_GEMs/stm_model_test/minimal_mouse_media_results")

for name, model in models.copy().items():
    
    min_media = dict(mouse_min_medium)
    tet_min_media = dict(mouse_min_medium)
    tet_min_media.update({
        "EX_tet_e": 4.0
    })

    min_model = model.copy()
    min_model.medium = min_media

    tet_model = model.copy()
    tet_model.medium = tet_min_media

    min_solution = min_model.optimize()
    tet_solution = tet_model.optimize()

    # Minimal media flux dataframe
    min_solution_flux_df = min_solution.fluxes.to_frame(name="Flux")

    min_solution_flux_df["Reaction Name"] = min_solution_flux_df.index.map(
        lambda rxn: min_model.reactions.get_by_id(rxn).name
    )

    min_solution_flux_df.index.name = "Reaction"

    min_solution_flux_df.to_csv(
        f"{results_dir}/{name}_minimal_media_solution_flux.csv"
    )

    # Tetrathionate minimal media flux dataframe
    tet_solution_flux_df = tet_solution.fluxes.to_frame(name="Flux")

    tet_solution_flux_df["Reaction Name"] = tet_solution_flux_df.index.map(
        lambda rxn: tet_model.reactions.get_by_id(rxn).name
    )

    tet_solution_flux_df.index.name = "Reaction"

    tet_solution_flux_df.to_csv(
        f"{results_dir}/{name}_tet_minimal_media_solution_flux.csv"
    )
    # Minimal media summary flux dataframe
    min_summary = min_model.summary(solution=min_solution)

    min_uptake = min_summary.uptake_flux.copy()
    min_uptake["Type"] = "Uptake"

    min_secretion = min_summary.secretion_flux.copy()
    min_secretion["Type"] = "Secretion"

    min_summary_flux_df = pd.concat(
        [min_uptake, min_secretion],
        ignore_index=True
    )   

    min_summary_flux_df["Objective Value"] = pd.NA
    min_summary_flux_df.loc[0, "Objective Value"] = min_solution.objective_value

    min_summary_flux_df.to_csv(
     f"{results_dir}/{name}_minimal_media_summary_flux.csv",
        index=False
    )

    # Tetrathionate minimal media summary flux dataframe
    tet_summary = tet_model.summary(solution=tet_solution)

    tet_uptake = tet_summary.uptake_flux.copy()
    tet_uptake["Type"] = "Uptake"

    tet_secretion = tet_summary.secretion_flux.copy()
    tet_secretion["Type"] = "Secretion"

    tet_summary_flux_df = pd.concat(
        [tet_uptake, tet_secretion],
        ignore_index=True
    )

    tet_summary_flux_df["Objective Value"] = pd.NA
    tet_summary_flux_df.loc[0, "Objective Value"] = tet_solution.objective_value

    tet_summary_flux_df.to_csv(
        f"{results_dir}/{name}_tet_minimal_media_summary_flux.csv",
        index=False
    )
    # Minimal media ATP summary dataframe
    min_atp_summary = min_model.metabolites.atp_c.summary(
        solution=min_solution
    )

    min_atp_producing = min_atp_summary.producing_flux.copy()
    min_atp_producing["Type"] = "Producing"

    min_atp_consuming = min_atp_summary.consuming_flux.copy()
    min_atp_consuming["Type"] = "Consuming"

    min_atp_summary_df = pd.concat(
        [min_atp_producing, min_atp_consuming],
        ignore_index=True
    )

    min_atp_summary_df.to_csv(
        f"{results_dir}/{name}_minimal_media_atp_summary.csv",
        index=False
    )


    # Tetrathionate minimal media ATP summary dataframe

    tet_atp_summary = tet_model.metabolites.atp_c.summary(
        solution=tet_solution
    )

    tet_atp_producing = tet_atp_summary.producing_flux.copy()
    tet_atp_producing["Type"] = "Producing"

    tet_atp_consuming = tet_atp_summary.consuming_flux.copy()
    tet_atp_consuming["Type"] = "Consuming"

    tet_atp_summary_df = pd.concat(
        [tet_atp_producing, tet_atp_consuming],
        ignore_index=True
    )

    tet_atp_summary_df.to_csv(
        f"{results_dir}/{name}_tet_minimal_media_atp_summary.csv",
        index=False
    )

    # Minimal media NADH summary dataframe

    min_nadh_summary = min_model.metabolites.nadh_c.summary(
        solution=min_solution
    )   

    min_nadh_producing = min_nadh_summary.producing_flux.copy()
    min_nadh_producing["Type"] = "Producing"

    min_nadh_consuming = min_nadh_summary.consuming_flux.copy()
    min_nadh_consuming["Type"] = "Consuming"

    min_nadh_summary_df = pd.concat(
        [min_nadh_producing, min_nadh_consuming],
        ignore_index=True
    )

    min_nadh_summary_df.to_csv(
        f"{results_dir}/{name}_minimal_media_nadh_summary.csv",
        index=False
    )

# Tetrathionate minimal media NADH summary dataframe

    tet_nadh_summary = tet_model.metabolites.nadh_c.summary(
        solution=tet_solution
    )

    tet_nadh_producing = tet_nadh_summary.producing_flux.copy()
    tet_nadh_producing["Type"] = "Producing"

    tet_nadh_consuming = tet_nadh_summary.consuming_flux.copy()
    tet_nadh_consuming["Type"] = "Consuming"

    tet_nadh_summary_df = pd.concat(
        [tet_nadh_producing, tet_nadh_consuming],
        ignore_index=True
    )

    tet_nadh_summary_df.to_csv(
        f"{results_dir}/{name}_tet_minimal_media_nadh_summary.csv",
        index=False
    )

for name, model in models.copy().items():

    min_media = dict(mouse_min_medium)

    succ_min_media = dict(mouse_min_medium)

    succ_min_media.update({
        "EX_succ_e": 3.0
    })

    min_model = model.copy()
    min_model.medium = min_media

    succ_model = model.copy()
    succ_model.medium = succ_min_media

    min_solution = min_model.optimize()
    succ_solution = succ_model.optimize()

    # Minimal media flux dataframe

    min_solution_flux_df = min_solution.fluxes.to_frame(name="Flux")

    min_solution_flux_df["Reaction Name"] = min_solution_flux_df.index.map(
        lambda rxn: min_model.reactions.get_by_id(rxn).name
    )

    min_solution_flux_df.index.name = "Reaction"

    min_solution_flux_df.to_csv(
        f"{results_dir}/{name}_minimal_media_solution_flux.csv"
    )

    # Succinate minimal media flux dataframe

    succ_solution_flux_df = succ_solution.fluxes.to_frame(name="Flux")

    succ_solution_flux_df["Reaction Name"] = succ_solution_flux_df.index.map(
        lambda rxn: succ_model.reactions.get_by_id(rxn).name
    )

    succ_solution_flux_df.index.name = "Reaction"

    succ_solution_flux_df.to_csv(
        f"{results_dir}/{name}_succ_minimal_media_solution_flux.csv"
    )

    # Minimal media summary flux dataframe

    min_summary = min_model.summary(solution=min_solution)

    min_uptake = min_summary.uptake_flux.copy()
    min_uptake["Type"] = "Uptake"

    min_secretion = min_summary.secretion_flux.copy()
    min_secretion["Type"] = "Secretion"

    min_summary_flux_df = pd.concat(
        [min_uptake, min_secretion],
        ignore_index=True
    )

    min_summary_flux_df["Objective Value"] = pd.NA
    min_summary_flux_df.loc[0, "Objective Value"] = min_solution.objective_value

    min_summary_flux_df.to_csv(
        f"{results_dir}/{name}_minimal_media_summary_flux.csv",
        index=False
    )

    # Succinate minimal media summary flux dataframe

    succ_summary = succ_model.summary(solution=succ_solution)

    succ_uptake = succ_summary.uptake_flux.copy()
    succ_uptake["Type"] = "Uptake"

    succ_secretion = succ_summary.secretion_flux.copy()
    succ_secretion["Type"] = "Secretion"

    succ_summary_flux_df = pd.concat(
        [succ_uptake, succ_secretion],
        ignore_index=True
    )

    succ_summary_flux_df["Objective Value"] = pd.NA
    succ_summary_flux_df.loc[0, "Objective Value"] = succ_solution.objective_value

    succ_summary_flux_df.to_csv(
        f"{results_dir}/{name}_succ_minimal_media_summary_flux.csv",
        index=False
    )

    # Minimal media ATP summary dataframe

    min_atp_summary = min_model.metabolites.atp_c.summary(
        solution=min_solution
    )

    min_atp_producing = min_atp_summary.producing_flux.copy()
    min_atp_producing["Type"] = "Producing"

    min_atp_consuming = min_atp_summary.consuming_flux.copy()
    min_atp_consuming["Type"] = "Consuming"

    min_atp_summary_df = pd.concat(
        [min_atp_producing, min_atp_consuming],
        ignore_index=True
    )

    min_atp_summary_df.to_csv(
        f"{results_dir}/{name}_minimal_media_atp_summary.csv",
        index=False
    )

    # Succinate minimal media ATP summary dataframe

    succ_atp_summary = succ_model.metabolites.atp_c.summary(
        solution=succ_solution
    )

    succ_atp_producing = succ_atp_summary.producing_flux.copy()
    succ_atp_producing["Type"] = "Producing"

    succ_atp_consuming = succ_atp_summary.consuming_flux.copy()
    succ_atp_consuming["Type"] = "Consuming"

    succ_atp_summary_df = pd.concat(
        [succ_atp_producing, succ_atp_consuming],
        ignore_index=True
    )

    succ_atp_summary_df.to_csv(
        f"{results_dir}/{name}_succ_minimal_media_atp_summary.csv",
        index=False
    )

    # Minimal media NADH summary dataframe

    min_nadh_summary = min_model.metabolites.nadh_c.summary(
        solution=min_solution
    )

    min_nadh_producing = min_nadh_summary.producing_flux.copy()
    min_nadh_producing["Type"] = "Producing"

    min_nadh_consuming = min_nadh_summary.consuming_flux.copy()
    min_nadh_consuming["Type"] = "Consuming"

    min_nadh_summary_df = pd.concat(
        [min_nadh_producing, min_nadh_consuming],
        ignore_index=True
    )

    min_nadh_summary_df.to_csv(
        f"{results_dir}/{name}_minimal_media_nadh_summary.csv",
        index=False
    )

    # Succinate minimal media NADH summary dataframe

    succ_nadh_summary = succ_model.metabolites.nadh_c.summary(
        solution=succ_solution
    )

    succ_nadh_producing = succ_nadh_summary.producing_flux.copy()
    succ_nadh_producing["Type"] = "Producing"

    succ_nadh_consuming = succ_nadh_summary.consuming_flux.copy()
    succ_nadh_consuming["Type"] = "Consuming"

    succ_nadh_summary_df = pd.concat(
        [succ_nadh_producing, succ_nadh_consuming],
        ignore_index=True
    )

    succ_nadh_summary_df.to_csv(
        f"{results_dir}/{name}_succ_minimal_media_nadh_summary.csv",
        index=False
    )

for name, model in models.copy().items():

    min_media = dict(mouse_min_medium)

    tet_and_succ_min_media = dict(mouse_min_medium)

    tet_and_succ_min_media.update({
        "EX_tet_e": 4.0,
        "EX_succ_e": 3.0
    })

    min_model = model.copy()
    min_model.medium = min_media

    tet_and_succ_model = model.copy()
    tet_and_succ_model.medium = tet_and_succ_min_media

    min_solution = min_model.optimize()
    tet_and_succ_solution = tet_and_succ_model.optimize()

    # Minimal media flux dataframe

    min_solution_flux_df = min_solution.fluxes.to_frame(name="Flux")

    min_solution_flux_df["Reaction Name"] = min_solution_flux_df.index.map(
        lambda rxn: min_model.reactions.get_by_id(rxn).name
    )

    min_solution_flux_df.index.name = "Reaction"

    min_solution_flux_df.to_csv(
        f"{results_dir}/{name}_minimal_media_solution_flux.csv"
    )

    # Tetrathionate and succinate minimal media flux dataframe

    tet_and_succ_solution_flux_df = tet_and_succ_solution.fluxes.to_frame(
        name="Flux"
    )

    tet_and_succ_solution_flux_df["Reaction Name"] = (
        tet_and_succ_solution_flux_df.index.map(
            lambda rxn: tet_and_succ_model.reactions.get_by_id(rxn).name
        )
    )

    tet_and_succ_solution_flux_df.index.name = "Reaction"

    tet_and_succ_solution_flux_df.to_csv(
        f"{results_dir}/{name}_tet_and_succ_minimal_media_solution_flux.csv"
    )

    # Minimal media summary flux dataframe

    min_summary = min_model.summary(solution=min_solution)

    min_uptake = min_summary.uptake_flux.copy()
    min_uptake["Type"] = "Uptake"

    min_secretion = min_summary.secretion_flux.copy()
    min_secretion["Type"] = "Secretion"

    min_summary_flux_df = pd.concat(
        [min_uptake, min_secretion],
        ignore_index=True
    )

    min_summary_flux_df["Objective Value"] = pd.NA
    min_summary_flux_df.loc[0, "Objective Value"] = min_solution.objective_value

    min_summary_flux_df.to_csv(
        f"{results_dir}/{name}_minimal_media_summary_flux.csv",
        index=False
    )

    # Tetrathionate and succinate minimal media summary flux dataframe

    tet_and_succ_summary = tet_and_succ_model.summary(
        solution=tet_and_succ_solution
    )

    tet_and_succ_uptake = tet_and_succ_summary.uptake_flux.copy()
    tet_and_succ_uptake["Type"] = "Uptake"

    tet_and_succ_secretion = tet_and_succ_summary.secretion_flux.copy()
    tet_and_succ_secretion["Type"] = "Secretion"

    tet_and_succ_summary_flux_df = pd.concat(
        [tet_and_succ_uptake, tet_and_succ_secretion],
        ignore_index=True
    )

    tet_and_succ_summary_flux_df["Objective Value"] = pd.NA
    tet_and_succ_summary_flux_df.loc[0, "Objective Value"] = (
        tet_and_succ_solution.objective_value
    )

    tet_and_succ_summary_flux_df.to_csv(
        f"{results_dir}/{name}_tet_and_succ_minimal_media_summary_flux.csv",
        index=False
    )

    # Minimal media ATP summary dataframe

    min_atp_summary = min_model.metabolites.atp_c.summary(
        solution=min_solution
    )

    min_atp_producing = min_atp_summary.producing_flux.copy()
    min_atp_producing["Type"] = "Producing"

    min_atp_consuming = min_atp_summary.consuming_flux.copy()
    min_atp_consuming["Type"] = "Consuming"

    min_atp_summary_df = pd.concat(
        [min_atp_producing, min_atp_consuming],
        ignore_index=True
    )

    min_atp_summary_df.to_csv(
        f"{results_dir}/{name}_minimal_media_atp_summary.csv",
        index=False
    )

    # Tetrathionate and succinate minimal media ATP summary dataframe

    tet_and_succ_atp_summary = tet_and_succ_model.metabolites.atp_c.summary(
        solution=tet_and_succ_solution
    )

    tet_and_succ_atp_producing = tet_and_succ_atp_summary.producing_flux.copy()
    tet_and_succ_atp_producing["Type"] = "Producing"

    tet_and_succ_atp_consuming = tet_and_succ_atp_summary.consuming_flux.copy()
    tet_and_succ_atp_consuming["Type"] = "Consuming"

    tet_and_succ_atp_summary_df = pd.concat(
        [tet_and_succ_atp_producing, tet_and_succ_atp_consuming],
        ignore_index=True
    )

    tet_and_succ_atp_summary_df.to_csv(
        f"{results_dir}/{name}_tet_and_succ_minimal_media_atp_summary.csv",
        index=False
    )

    # Minimal media NADH summary dataframe

    min_nadh_summary = min_model.metabolites.nadh_c.summary(
        solution=min_solution
    )

    min_nadh_producing = min_nadh_summary.producing_flux.copy()
    min_nadh_producing["Type"] = "Producing"

    min_nadh_consuming = min_nadh_summary.consuming_flux.copy()
    min_nadh_consuming["Type"] = "Consuming"

    min_nadh_summary_df = pd.concat(
        [min_nadh_producing, min_nadh_consuming],
        ignore_index=True
    )

    min_nadh_summary_df.to_csv(
        f"{results_dir}/{name}_minimal_media_nadh_summary.csv",
        index=False
    )

    # Tetrathionate and succinate minimal media NADH summary dataframe

    tet_and_succ_nadh_summary = tet_and_succ_model.metabolites.nadh_c.summary(
        solution=tet_and_succ_solution
    )

    tet_and_succ_nadh_producing = (
        tet_and_succ_nadh_summary.producing_flux.copy()
    )
    tet_and_succ_nadh_producing["Type"] = "Producing"

    tet_and_succ_nadh_consuming = (
        tet_and_succ_nadh_summary.consuming_flux.copy()
    )
    tet_and_succ_nadh_consuming["Type"] = "Consuming"

    tet_and_succ_nadh_summary_df = pd.concat(
        [tet_and_succ_nadh_producing, tet_and_succ_nadh_consuming],
        ignore_index=True
    )

    tet_and_succ_nadh_summary_df.to_csv(
        f"{results_dir}/{name}_tet_and_succ_minimal_media_nadh_summary.csv",
        index=False
    )

for name, model in models.copy().items():

    min_media = dict(mouse_min_medium)

    ppd_min_media = dict(mouse_min_medium)

    ppd_min_media.update({
        "EX_12ppd__R_e": 0.1,
        "EX_12ppd__S_e": 0.1
    })

    min_model = model.copy()
    min_model.medium = min_media

    ppd_model = model.copy()
    ppd_model.medium = ppd_min_media

    min_solution = min_model.optimize()
    ppd_solution = ppd_model.optimize()

    # Minimal media flux dataframe

    min_solution_flux_df = min_solution.fluxes.to_frame(name="Flux")

    min_solution_flux_df["Reaction Name"] = min_solution_flux_df.index.map(
        lambda rxn: min_model.reactions.get_by_id(rxn).name
    )

    min_solution_flux_df.index.name = "Reaction"

    min_solution_flux_df.to_csv(
        f"{results_dir}/{name}_minimal_media_solution_flux.csv"
    )

    # PPD minimal media flux dataframe

    ppd_solution_flux_df = ppd_solution.fluxes.to_frame(name="Flux")

    ppd_solution_flux_df["Reaction Name"] = ppd_solution_flux_df.index.map(
        lambda rxn: ppd_model.reactions.get_by_id(rxn).name
    )

    ppd_solution_flux_df.index.name = "Reaction"

    ppd_solution_flux_df.to_csv(
        f"{results_dir}/{name}_ppd_minimal_media_solution_flux.csv"
    )

    # Minimal media summary flux dataframe

    min_summary = min_model.summary(solution=min_solution)

    min_uptake = min_summary.uptake_flux.copy()
    min_uptake["Type"] = "Uptake"

    min_secretion = min_summary.secretion_flux.copy()
    min_secretion["Type"] = "Secretion"

    min_summary_flux_df = pd.concat(
        [min_uptake, min_secretion],
        ignore_index=True
    )

    min_summary_flux_df["Objective Value"] = pd.NA
    min_summary_flux_df.loc[0, "Objective Value"] = min_solution.objective_value

    min_summary_flux_df.to_csv(
        f"{results_dir}/{name}_minimal_media_summary_flux.csv",
        index=False
    )

    # PPD minimal media summary flux dataframe

    ppd_summary = ppd_model.summary(solution=ppd_solution)

    ppd_uptake = ppd_summary.uptake_flux.copy()
    ppd_uptake["Type"] = "Uptake"

    ppd_secretion = ppd_summary.secretion_flux.copy()
    ppd_secretion["Type"] = "Secretion"

    ppd_summary_flux_df = pd.concat(
        [ppd_uptake, ppd_secretion],
        ignore_index=True
    )

    ppd_summary_flux_df["Objective Value"] = pd.NA
    ppd_summary_flux_df.loc[0, "Objective Value"] = ppd_solution.objective_value

    ppd_summary_flux_df.to_csv(
        f"{results_dir}/{name}_ppd_minimal_media_summary_flux.csv",
        index=False
    )

    # Minimal media ATP summary dataframe

    min_atp_summary = min_model.metabolites.atp_c.summary(
        solution=min_solution
    )

    min_atp_producing = min_atp_summary.producing_flux.copy()
    min_atp_producing["Type"] = "Producing"

    min_atp_consuming = min_atp_summary.consuming_flux.copy()
    min_atp_consuming["Type"] = "Consuming"

    min_atp_summary_df = pd.concat(
        [min_atp_producing, min_atp_consuming],
        ignore_index=True
    )

    min_atp_summary_df.to_csv(
        f"{results_dir}/{name}_minimal_media_atp_summary.csv",
        index=False
    )

    # PPD minimal media ATP summary dataframe

    ppd_atp_summary = ppd_model.metabolites.atp_c.summary(
        solution=ppd_solution
    )

    ppd_atp_producing = ppd_atp_summary.producing_flux.copy()
    ppd_atp_producing["Type"] = "Producing"

    ppd_atp_consuming = ppd_atp_summary.consuming_flux.copy()
    ppd_atp_consuming["Type"] = "Consuming"

    ppd_atp_summary_df = pd.concat(
        [ppd_atp_producing, ppd_atp_consuming],
        ignore_index=True
    )

    ppd_atp_summary_df.to_csv(
        f"{results_dir}/{name}_ppd_minimal_media_atp_summary.csv",
        index=False
    )

    # Minimal media NADH summary dataframe

    min_nadh_summary = min_model.metabolites.nadh_c.summary(
        solution=min_solution
    )

    min_nadh_producing = min_nadh_summary.producing_flux.copy()
    min_nadh_producing["Type"] = "Producing"

    min_nadh_consuming = min_nadh_summary.consuming_flux.copy()
    min_nadh_consuming["Type"] = "Consuming"

    min_nadh_summary_df = pd.concat(
        [min_nadh_producing, min_nadh_consuming],
        ignore_index=True
    )

    min_nadh_summary_df.to_csv(
        f"{results_dir}/{name}_minimal_media_nadh_summary.csv",
        index=False
    )

    # PPD minimal media NADH summary dataframe

    ppd_nadh_summary = ppd_model.metabolites.nadh_c.summary(
        solution=ppd_solution
    )

    ppd_nadh_producing = ppd_nadh_summary.producing_flux.copy()
    ppd_nadh_producing["Type"] = "Producing"

    ppd_nadh_consuming = ppd_nadh_summary.consuming_flux.copy()
    ppd_nadh_consuming["Type"] = "Consuming"

    ppd_nadh_summary_df = pd.concat(
        [ppd_nadh_producing, ppd_nadh_consuming],
        ignore_index=True
    )

    ppd_nadh_summary_df.to_csv(
        f"{results_dir}/{name}_ppd_minimal_media_nadh_summary.csv",
        index=False
    )

for name, model in models.copy().items():

    tet_succ_ppd_min_media = dict(mouse_min_medium)

    tet_succ_ppd_min_media.update({
        "EX_succ_e": 3.0,
        "EX_tet_e": 4.0,
        "EX_12ppd__R_e": 0.1,
        "EX_12ppd__S_e": 0.1
    })

    tet_succ_ppd_model = model.copy()

    tet_succ_ppd_model.medium = tet_succ_ppd_min_media

    tet_succ_ppd_solution = tet_succ_ppd_model.optimize()

# Tetrathionate, succinate, and PPD minimal media flux dataframe

    tet_succ_ppd_solution_flux_df = tet_succ_ppd_solution.fluxes.to_frame(
        name="Flux"
    )

    tet_succ_ppd_solution_flux_df["Reaction Name"] = (
        tet_succ_ppd_solution_flux_df.index.map(
            lambda rxn: tet_succ_ppd_model.reactions.get_by_id(rxn).name
        )
    )

    tet_succ_ppd_solution_flux_df.index.name = "Reaction"

    tet_succ_ppd_solution_flux_df.to_csv(
        f"{results_dir}/{name}_tet_succ_ppd_minimal_media_solution_flux.csv"
    )

# Tetrathionate, succinate, and PPD minimal media summary flux dataframe

    tet_succ_ppd_summary = tet_succ_ppd_model.summary(
        solution=tet_succ_ppd_solution
    )

    tet_succ_ppd_uptake = tet_succ_ppd_summary.uptake_flux.copy()

    tet_succ_ppd_uptake["Type"] = "Uptake"

    tet_succ_ppd_secretion = tet_succ_ppd_summary.secretion_flux.copy()

    tet_succ_ppd_secretion["Type"] = "Secretion"

    tet_succ_ppd_summary_flux_df = pd.concat(
        [tet_succ_ppd_uptake, tet_succ_ppd_secretion],
        ignore_index=True
    )

    tet_succ_ppd_summary_flux_df["Objective Value"] = pd.NA

    tet_succ_ppd_summary_flux_df.loc[0, "Objective Value"] = (
        tet_succ_ppd_solution.objective_value
    )

    tet_succ_ppd_summary_flux_df.to_csv(
        f"{results_dir}/{name}_tet_succ_ppd_minimal_media_summary_flux.csv",
        index=False
    )

# Tetrathionate, succinate, and PPD minimal media ATP summary dataframe

    tet_succ_ppd_atp_summary = tet_succ_ppd_model.metabolites.atp_c.summary(
        solution=tet_succ_ppd_solution
    )

    tet_succ_ppd_atp_producing = (
        tet_succ_ppd_atp_summary.producing_flux.copy()
    )

    tet_succ_ppd_atp_producing["Type"] = "Producing"

    tet_succ_ppd_atp_consuming = (
        tet_succ_ppd_atp_summary.consuming_flux.copy()
    )

    tet_succ_ppd_atp_consuming["Type"] = "Consuming"

    tet_succ_ppd_atp_summary_df = pd.concat(
        [tet_succ_ppd_atp_producing, tet_succ_ppd_atp_consuming],
        ignore_index=True
    )

    tet_succ_ppd_atp_summary_df.to_csv(
        f"{results_dir}/{name}_tet_succ_ppd_minimal_media_atp_summary.csv",
        index=False
    )

# Tetrathionate, succinate, and PPD minimal media NADH summary dataframe

    tet_succ_ppd_nadh_summary = tet_succ_ppd_model.metabolites.nadh_c.summary(
        solution=tet_succ_ppd_solution
    )

    tet_succ_ppd_nadh_producing = (
        tet_succ_ppd_nadh_summary.producing_flux.copy()
    )

    tet_succ_ppd_nadh_producing["Type"] = "Producing"

    tet_succ_ppd_nadh_consuming = (
        tet_succ_ppd_nadh_summary.consuming_flux.copy()
    )

    tet_succ_ppd_nadh_consuming["Type"] = "Consuming"

    tet_succ_ppd_nadh_summary_df = pd.concat(
        [tet_succ_ppd_nadh_producing, tet_succ_ppd_nadh_consuming],
        ignore_index=True
    )

    tet_succ_ppd_nadh_summary_df.to_csv(
        f"{results_dir}/{name}_tet_succ_ppd_minimal_media_nadh_summary.csv",
        index=False
    )