
# Set the gurobi solver 

import os
from pathlib import Path
from cobra.core.configuration import Configuration

# Use an existing cluster license if present; otherwise use the local Mac license file if it exists.
local_gurobi_license = Path.home() / ".gurobi.lic"
if local_gurobi_license.exists():
    os.environ.setdefault("GRB_LICENSE_FILE", str(local_gurobi_license))

Configuration().solver = "gurobi"


# In[5]:


# Read in the gene weights CSV files

import pandas as pd

stm_threshold_10_90 = pd.read_csv("../imat_prep_test/stm_threshold_10_90.csv")
stm_threshold_15_85 = pd.read_csv("../imat_prep_test/stm_threshold_15_85.csv")
stm_threshold_25_75 = pd.read_csv("../imat_prep_test/stm_threshold_25_75.csv")

mixed_stm_threshold_10_90 = pd.read_csv("../imat_prep_test/mixed_stm_threshold_10_90.csv")
mixed_stm_threshold_15_85 = pd.read_csv("../imat_prep_test/mixed_stm_threshold_15_85.csv")
mixed_stm_threshold_25_75 = pd.read_csv("../imat_prep_test/mixed_stm_threshold_25_75.csv")


# In[6]:


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


# In[9]:


# Load STM_V1_0

from cobra.io import load_model, load_json_model

control_model = load_model("STM_v1_0")


# In[10]:


# Save all genes in the control model in a set

model_genes = set([g.id for g in control_model.genes])

# QC: Check the count of genes
print("STM_V1_0 Control Model has", len(model_genes), "total genes")


# In[11]:


# Subset the gene weights lists for each sample to keep only the genes also found in the control model

stm1090_model = stm_threshold_10_90[stm_threshold_10_90["genes"].isin(model_genes)]
stm1585_model = stm_threshold_15_85[stm_threshold_15_85["genes"].isin(model_genes)]
stm2575_model = stm_threshold_25_75[stm_threshold_25_75["genes"].isin(model_genes)]

mixed_stm1090_model = mixed_stm_threshold_10_90[mixed_stm_threshold_10_90["genes"].isin(model_genes)]
mixed_stm1585_model = mixed_stm_threshold_15_85[mixed_stm_threshold_15_85["genes"].isin(model_genes)]
mixed_stm2575_model = mixed_stm_threshold_25_75[mixed_stm_threshold_25_75["genes"].isin(model_genes)]


# In[12]:


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
    subset_weights_mixedstm1585,
    subset_weights_mixedstm2575]

# 2. Make a list of how many genes are in each list
gene_counts = [len(weights) for weights in subset_gene_weights]

# 3. Compare how many genes are in each list and print the total if they are all the same
if len(set(gene_counts)) == 1:
    print(f"All models have {gene_counts[0]} genes")
else:
    print("Gene counts differ between models")
    print(gene_counts)


# ### Convert Gene Weights to Reaction Weights

# In[14]:


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


# ### Generate Constrained Models with iMAT

# In[15]:


# Generate constrained models with iMAT by integrating the reaction weights

from imatpy.model_creation import generate_model

stm1090_model = generate_model(model = control_model, rxn_weights = rxn_weights_stm1090, method = "imat_restrictions")
stm1585_model = generate_model(model = control_model, rxn_weights = rxn_weights_stm1585, method = "imat_restrictions")
stm2575_model = generate_model(model = control_model, rxn_weights = rxn_weights_stm2575, method = "imat_restrictions")

mixed_stm1090_model = generate_model(model = control_model, rxn_weights = rxn_weights_mixedstm1090, method = "imat_restrictions")
mixed_stm1585_model = generate_model(model = control_model, rxn_weights = rxn_weights_mixedstm1585, method = "imat_restrictions")
mixed_stm2575_model = generate_model(model = control_model, rxn_weights = rxn_weights_mixedstm2575, method = "imat_restrictions")


# In[16]:


# List the default medium exchange reactions, metabolites, and uptake value

for rxn_id, uptake_value in control_model.medium.items():
    rxn = control_model.reactions.get_by_id(rxn_id)
    met = list(rxn.metabolites)[0]
    print(rxn_id, "|", met.id, "|", met.name, "|", uptake_value)


# In[17]:


# Create a default_medium dictionary

default_medium = {
 'EX_ca2_e': 1000.0,
 'EX_cbl1_e': 0.01,
 'EX_cl_e': 1000.0,
 'EX_co2_e': 1000.0,
 'EX_cobalt2_e': 1000.0,
 'EX_cu2_e': 1000.0,
 'EX_mg2_e': 1000.0,
 'EX_fe2_e': 1000.0,
 'EX_fe3_e': 1000.0,
 'EX_mn2_e': 1000.0,
 'EX_mobd_e': 1000.0,
 'EX_na1_e': 1000.0,
 'EX_nh4_e': 1000.0,
 'EX_o2_e': 18.5,
 'EX_pi_e': 1000.0,
 'EX_zn2_e': 1000.0,
 'EX_so4_e': 1000.0,
 'EX_glc__D_e': 5.0,
 'EX_tungs_e': 1000.0,
 'EX_h2o_e': 100.0,
 'EX_h_e': 100.0,
 'EX_k_e': 1000.0}


# In[18]:


# Make an infection condition medium dictionary with metabolites from SW's list and include the default medium conditions.

infected_medium = dict(default_medium)
infected_medium.update({
    "EX_o2_e": 9.0,
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
    "EX_glc__D_e": 0.0, # eliminate the glucose present in the default medium
})


# In[19]:


# Make a non-infection condition medium dictionary with metabolites from SW's list and include the default medium conditions.

uninfected_medium = dict(default_medium)
uninfected_medium.update({
    "EX_o2_e": 9.0,
    "EX_tet_e": 0.05,
    "EX_no3_e": 0.05,
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
    "EX_glc__D_e": 0.0, # eliminate the glucose present in the default medium
})


# ## Oxygen Gradient Study
# 
# ### Making a function that runs the simulation
# 
# 1. Make a dictionary with the medium conditions and names as the keys.
# 
# I need a function that has arguments for the model name, medium conditon, and increment/step value to increase oxygen level.
# 
# Here's what I need it to do:
# 
# - Take the `step_value` argument and make a list of the `oxygen_levels` which are the number increment values of oxygen flux that make up the gradient.
# - Take the medium condition, copy, and make a new copy for each oxygen lecel in the `oxygen_levels` list.
# - Take the model name, copy the model, and make a new copy with each new medium condition for each oxygen level.
# 
# Now, I should have X amount of copies of the model equal to the total amount of oxygen increments in the oxygen_levels list. 
# 
# - Run FBA on all the models.
# - Put the results in a dictionary.
# - Print as a dataframe. 

# In[22]:


base_media = {
    "default": default_medium,
    "infected": infected_medium,
    "uninfected": uninfected_medium,
}


# In[25]:


original_models = {
    "control": control_model,
    "stm1090": stm1090_model,
    "stm1585": stm1585_model,
    "stm2575": stm2575_model,
    "mixed_stm1090": mixed_stm1090_model,
    "mixed_stm1585": mixed_stm1585_model,
    "mixed_stm2575": mixed_stm2575_model,
}


# In[26]:


fba_models = {
    model_name: model.copy()
    for model_name, model in original_models.items()
}


# In[21]:


import numpy as np
import pandas as pd

def run_oxygen_sweep(model_name, medium_condition, step_level=0.5,
                     o2_min=0.0, o2_max=18.5, o2_rxn="EX_o2_e"):
    """
    Run FBA across an adjustable oxygen gradient for one model and one medium condition.

    Parameters
    ----------
    model_name : str
        Name of the model to use, such as "control", "stm1090", or "mixed_stm2575".

    medium_condition : str
        Name of the base medium condition, such as "default", "infected", or "uninfected".

    step_level : float
        Oxygen increment size. Examples: 0.1, 0.5, 1, or 5.

    o2_min : float
        Starting oxygen uptake value. Default is 0.0.

    o2_max : float
        Maximum oxygen uptake value. Default is 18.5.

    o2_rxn : str
        Exchange reaction ID for oxygen. Default is "EX_o2_e".

    Returns
    -------
    pandas.DataFrame
        DataFrame with model name, medium condition, oxygen level, solver status,
        and objective value.
    """

    # Check that the requested medium exists
    if medium_condition not in base_media:
        raise ValueError(
            f"{medium_condition} is not in base_media. "
            f"Available media: {list(base_media.keys())}"
        )

    # Check that the requested model exists
    if model_name not in fba_models:
        raise ValueError(
            f"{model_name} is not in fba_models. "
            f"Available models: {list(fba_models.keys())}"
        )

    # Check that the oxygen step is valid
    if step_level <= 0:
        raise ValueError("step_level must be greater than 0.")

    # Build oxygen levels from o2_min to o2_max using the selected step size
    oxygen_levels = np.arange(o2_min, o2_max + step_level, step_level)

    # Remove values that go above o2_max
    oxygen_levels = oxygen_levels[oxygen_levels <= o2_max]

    # Round values to avoid floating-point artifacts like 0.30000000004
    oxygen_levels = np.round(oxygen_levels, 6)

    # Make sure o2_max is included even if the step does not land on it exactly
    if oxygen_levels[-1] != o2_max:
        oxygen_levels = np.append(oxygen_levels, o2_max)

    # Get the selected base medium
    base_medium = base_media[medium_condition]

    # Get a fresh copy of the selected model
    model = fba_models[model_name].copy()

    results = []

    # Loop through each oxygen value
    for o2_value in oxygen_levels:

        # Copy the base medium so the original dictionary is not changed
        medium = base_medium.copy()

        # Set oxygen to the current oxygen value
        medium[o2_rxn] = float(o2_value)

        try:
            with model:
                model.medium = medium
                solution = model.optimize()

                if solution.status == "optimal":
                    objective_value = solution.objective_value
                else:
                    objective_value = np.nan

                results.append({
                    "model": model_name,
                    "medium_condition": medium_condition,
                    "oxygen": float(o2_value),
                    "status": solution.status,
                    "objective_value": objective_value
                })

        except Exception as error:
            results.append({
                "model": model_name,
                "medium_condition": medium_condition,
                "oxygen": float(o2_value),
                "status": "failed",
                "objective_value": np.nan,
                "error": str(error)
            })

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("oxygen").reset_index(drop=True)

    return results_df


# ## Preparing Runs

# ## Control Run

# In[32]:


control_models = [
    "control"
]

control_results = []

for model_name in control_models:
    for medium_condition in base_media:

        print(f"Running {model_name} in {medium_condition} medium")

        try:
            results_df = run_oxygen_sweep(
                model_name=model_name,
                medium_condition=medium_condition,
                step_level=0.5
            )

            results_df["sample_type"] = "control"
            control_results.append(results_df)

        except Exception as error:
            print(f"FAILED: {model_name} in {medium_condition}")
            print(error)

            failed_df = pd.DataFrame([{
                "sample_type": "control",
                "model": model_name,
                "medium_condition": medium_condition,
                "oxygen": np.nan,
                "status": "failed_chunk",
                "objective_value": np.nan,
                "error": str(error)
            }])

            control_results.append(failed_df)
            continue

control_oxygen_sweep_df = pd.concat(
    control_results,
    ignore_index=True
)


# In[34]:


control_oxygen_sweep_df.head()


# ## STM-Only Run

# In[30]:


stm_models = [
    "stm1090",
    "stm1585",
    "stm2575"
]

stm_results = []

for model_name in stm_models:
    for medium_condition in base_media:

        print(f"Running {model_name} in {medium_condition} medium")

        try:
            results_df = run_oxygen_sweep(
                model_name=model_name,
                medium_condition=medium_condition,
                step_level=0.5
            )

            results_df["sample_type"] = "stm"
            stm_results.append(results_df)

        except Exception as error:
            print(f"FAILED: {model_name} in {medium_condition}")
            print(error)

            failed_df = pd.DataFrame([{
                "sample_type": "stm",
                "model": model_name,
                "medium_condition": medium_condition,
                "oxygen": np.nan,
                "status": "failed_chunk",
                "objective_value": np.nan,
                "error": str(error)
            }])

            stm_results.append(failed_df)
            continue

stm_oxygen_sweep_df = pd.concat(
    stm_results,
    ignore_index=True
)

stm_oxygen_sweep_df.head()


# ## Mixed STM Run

# In[37]:


print(type(mixed_stm1090_model))
print(type(mixed_stm1585_model))
print(type(mixed_stm2575_model))


# In[31]:


# Mixed STM oxygen gradient results

mixed_stm_models = [
    "mixed_stm1090",
    "mixed_stm1585",
    "mixed_stm2575"
]

mixed_stm_results = []

for model_name in mixed_stm_models:
    for medium_condition in base_media:

        print(f"Running {model_name} in {medium_condition} medium")

        try:
            results_df = run_oxygen_sweep(
                model_name=model_name,
                medium_condition=medium_condition,
                step_level=0.5
            )

            results_df["sample_type"] = "mixed_stm"
            mixed_stm_results.append(results_df)

        except Exception as error:
            print(f"FAILED: {model_name} in {medium_condition}")
            print(error)

            failed_df = pd.DataFrame([{
                "sample_type": "mixed_stm",
                "model": model_name,
                "medium_condition": medium_condition,
                "oxygen": np.nan,
                "status": "failed_chunk",
                "objective_value": np.nan,
                "error": str(error)
            }])

            mixed_stm_results.append(failed_df)
            continue

mixed_stm_oxygen_sweep_df = pd.concat(
    mixed_stm_results,
    ignore_index=True
)

mixed_stm_oxygen_sweep_df.head()


# ## Save Results!

# In[36]:


import os

os.makedirs("oxygen_sweep_results_062626", exist_ok=True)

control_oxygen_sweep_df.to_csv("oxygen_sweep_results_062626/control_oxygen_sweep_df.csv", index=False)
stm_oxygen_sweep_df.to_csv("oxygen_sweep_results_062626/stm_oxygen_sweep_df.csv", index=False)
mixed_stm_oxygen_sweep_df.to_csv("oxygen_sweep_results_062626/mixed_stm_oxygen_sweep_df.csv", index=False)


# ## Ploting
# 
# Combine all the result dataframes. Make a plotting function. plot all. 

# In[ ]:


# combine all the result dataframes 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Relative output directory
output_dir = Path("062226_labmeeting")

# Create the directory if it does not already exist
output_dir.mkdir(parents=True, exist_ok=True)

all_oxygen_sweep_df = pd.concat(
    [
        control_oxygen_sweep_df,
        stm_oxygen_sweep_df,
        mixed_stm_oxygen_sweep_df
    ],
    ignore_index=True
)


# In[38]:


def save_oxygen_sweep_plot(
    df,
    filter_column,
    filter_value,
    group_column,
    title,
    filename,
    xlim=(0, 18.5),
    ylim=(0, 1.8)
):
    """
    Plot and save oxygen sweep results.
    """

    plot_df = df[df[filter_column] == filter_value].copy()

    # Remove complete failed_chunk rows with no oxygen value
    plot_df = plot_df.dropna(subset=["oxygen"])

    # Treat failed or missing objective values as no growth
    plot_df["plot_objective_value"] = plot_df["objective_value"].fillna(0)

    plt.figure(figsize=(8, 5))

    for group_name, group_df in plot_df.groupby(group_column):

        group_df = group_df.sort_values("oxygen")

        plt.plot(
            group_df["oxygen"],
            group_df["plot_objective_value"],
            marker="o",
            linewidth=1,
            markersize=2,
            label=group_name
        )

        # Mark failed, missing, or zero-growth points in red
        no_growth_df = group_df[
            (group_df["objective_value"].isna()) |
            (np.isclose(group_df["plot_objective_value"], 0.0)) |
            (group_df["status"] != "optimal")
        ]

        plt.scatter(
            no_growth_df["oxygen"],
            no_growth_df["plot_objective_value"],
            color="red",
            s=40,
            zorder=3
        )

    plt.xlim(xlim)

    if ylim is not None:
        plt.ylim(ylim)

    plt.xlabel("Oxygen uptake value for EX_o2_e")
    plt.ylabel("Objective value")
    plt.title(title)
    plt.legend()
    plt.grid(True)

    # Save plot
    save_path = output_dir / filename
    plt.savefig(save_path, dpi=300, bbox_inches="tight")

    plt.show()

    print(f"Saved: {save_path}")


# ### Control model with three medium conditions

# In[ ]:


save_oxygen_sweep_plot(
    df=all_oxygen_sweep_df,
    filter_column="model",
    filter_value="control",
    group_column="medium_condition",
    title="Control model growth across oxygen gradient",
    filename="plot_01_control_three_media.png"
)


# ### stm1090 with three medium conditions

# In[ ]:


save_oxygen_sweep_plot(
    df=all_oxygen_sweep_df,
    filter_column="model",
    filter_value="stm1090",
    group_column="medium_condition",
    title="STM 10/90 model growth across oxygen gradient",
    filename="plot_02_stm1090_three_media.png"
)


# ### mixed_stm1090 with three medium conditions

# In[ ]:


save_oxygen_sweep_plot(
    df=all_oxygen_sweep_df,
    filter_column="model",
    filter_value="mixed_stm1090",
    group_column="medium_condition",
    title="Mixed STM 10/90 model growth across oxygen gradient",
    filename="plot_03_mixed_stm1090_three_media.png"
)


# ### Make a datafrmae of only the 10/90 models to compare control, stm_190, and mixed_stm1090 in each medium

# In[ ]:


selected_models = [
    "control",
    "stm1090",
    "mixed_stm1090"
]

selected_model_df = all_oxygen_sweep_df[
    all_oxygen_sweep_df["model"].isin(selected_models)
].copy()


# ### control, stm1090, and mixed_stm1090 in default medium

# In[ ]:


save_oxygen_sweep_plot(
    df=selected_model_df,
    filter_column="medium_condition",
    filter_value="default",
    group_column="model",
    title="Control, STM 10/90, and mixed STM 10/90 growth in default medium",
    filename="plot_04_three_models_default_medium.png"
)


# ### control, stm1090, and mixed_stm1090 in infected medium

# In[ ]:


save_oxygen_sweep_plot(
    df=selected_model_df,
    filter_column="medium_condition",
    filter_value="infected",
    group_column="model",
    title="Control, STM 10/90, and mixed STM 10/90 growth in infected medium",
    filename="plot_05_three_models_infected_medium.png"
)


# ### control, stm1090, and mixed_stm1090 in uninfected medium

# In[ ]:


save_oxygen_sweep_plot(
    df=selected_model_df,
    filter_column="medium_condition",
    filter_value="uninfected",
    group_column="model",
    title="Control, STM 10/90, and mixed STM 10/90 growth in uninfected medium",
    filename="plot_06_three_models_uninfected_medium.png"
)

