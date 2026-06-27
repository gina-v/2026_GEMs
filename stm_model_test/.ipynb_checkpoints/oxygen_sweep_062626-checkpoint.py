#!/usr/bin/env python
# coding: utf-8

"""
Oxygen sweep study for STM GEMs.

This script generates control, STM-only, and mixed-STM iMAT-constrained models,
runs FBA across an oxygen gradient, and saves all CSV and plot outputs to:

    <stm_model_test>/oxygen_sweep_results/062626

The script is written to run from either a local computer or the Farm cluster.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib

# Use a non-interactive backend so plots can be saved on the cluster/headless jobs.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cobra.core.configuration import Configuration
from cobra.io import load_model
from imatpy.model_creation import generate_model
from imatpy.parse_gpr import gene_to_rxn_weights


# -----------------------------------------------------------------------------
# Paths and solver setup
# -----------------------------------------------------------------------------

def find_stm_model_test_dir() -> Path:
    """Find the stm_model_test project directory on either local or cluster."""

    # Optional override for either local or cluster runs.
    env_path = os.environ.get("STM_MODEL_TEST_DIR")
    if env_path:
        return Path(env_path).expanduser().resolve()

    cwd = Path.cwd().resolve()

    # Use the current directory or one of its parents if it is stm_model_test.
    for path in [cwd, *cwd.parents]:
        if path.name == "stm_model_test":
            return path

    # Fall back to the expected project location under the user's home directory.
    fallback_path = Path.home() / "2026_GEMs" / "stm_model_test"
    if fallback_path.exists():
        return fallback_path.resolve()

    raise FileNotFoundError(
        "Could not find the stm_model_test directory. Run this script from inside "
        "stm_model_test, or set STM_MODEL_TEST_DIR to the full path."
    )


project_dir = find_stm_model_test_dir()
input_dir = project_dir.parent / "imat_prep_test"
output_dir = project_dir / "oxygen_sweep_results" / "062626"
output_dir.mkdir(parents=True, exist_ok=True)

print(f"Project directory: {project_dir}")
print(f"Input directory:   {input_dir}")
print(f"Output directory:  {output_dir}")

# Use an existing cluster license if present. Otherwise, use the local Mac license
# file if it exists. This avoids hard-coding a Mac-only path.
local_gurobi_license = Path.home() / ".gurobi.lic"
if local_gurobi_license.exists():
    os.environ.setdefault("GRB_LICENSE_FILE", str(local_gurobi_license))

# Keep Gurobi as the default solver, but allow override with COBRA_SOLVER.
Configuration().solver = os.environ.get("COBRA_SOLVER", "gurobi")
print(f"COBRA solver: {Configuration().solver}")


# -----------------------------------------------------------------------------
# Gene weights and model generation
# -----------------------------------------------------------------------------

def read_threshold_csv(filename: str) -> pd.DataFrame:
    """Read a threshold CSV from the imat_prep_test directory."""

    file_path = input_dir / filename
    if not file_path.exists():
        raise FileNotFoundError(f"Missing input file: {file_path}")

    return pd.read_csv(file_path)


# Read threshold CSV files for STM-only and mixed-STM samples.
stm_threshold_10_90 = read_threshold_csv("stm_threshold_10_90.csv")
stm_threshold_15_85 = read_threshold_csv("stm_threshold_15_85.csv")
stm_threshold_25_75 = read_threshold_csv("stm_threshold_25_75.csv")

mixed_stm_threshold_10_90 = read_threshold_csv("mixed_stm_threshold_10_90.csv")
mixed_stm_threshold_15_85 = read_threshold_csv("mixed_stm_threshold_15_85.csv")
mixed_stm_threshold_25_75 = read_threshold_csv("mixed_stm_threshold_25_75.csv")


# QC: Count genes in each threshold CSV before subsetting to model genes.
stm_set_10_90 = set(stm_threshold_10_90["genes"].astype(str).str.strip())
stm_set_15_85 = set(stm_threshold_15_85["genes"].astype(str).str.strip())
stm_set_25_75 = set(stm_threshold_25_75["genes"].astype(str).str.strip())

mixed_stm_set_10_90 = set(mixed_stm_threshold_10_90["genes"].astype(str).str.strip())
mixed_stm_set_15_85 = set(mixed_stm_threshold_15_85["genes"].astype(str).str.strip())
mixed_stm_set_25_75 = set(mixed_stm_threshold_25_75["genes"].astype(str).str.strip())

print(
    "STM gene counts:",
    f"10/90={len(stm_set_10_90)}",
    f"15/85={len(stm_set_15_85)}",
    f"25/75={len(stm_set_25_75)}",
)
print(
    "Mixed STM gene counts:",
    f"10/90={len(mixed_stm_set_10_90)}",
    f"15/85={len(mixed_stm_set_15_85)}",
    f"25/75={len(mixed_stm_set_25_75)}",
)


# Load the control model from COBRApy's model repository/cache.
control_model = load_model("STM_v1_0")
model_genes = sorted(g.id for g in control_model.genes)
model_gene_set = set(model_genes)
print(f"STM_V1_0 control model genes: {len(model_genes)}")


def align_gene_weights(threshold_df: pd.DataFrame, sample_name: str) -> pd.Series:
    """Keep model genes only, then align gene weights to all control-model genes."""

    required_columns = {"genes", "threshold_category"}
    missing_columns = required_columns - set(threshold_df.columns)
    if missing_columns:
        raise ValueError(f"{sample_name} is missing columns: {missing_columns}")

    threshold_df = threshold_df.copy()
    threshold_df["genes"] = threshold_df["genes"].astype(str).str.strip()

    model_gene_weights = threshold_df[threshold_df["genes"].isin(model_gene_set)]
    aligned_weights = (
        model_gene_weights.set_index("genes")["threshold_category"]
        .reindex(model_genes, fill_value=0)
    )

    return aligned_weights


# Align all sample gene weights to the same control-model gene list.
subset_weights_stm1090 = align_gene_weights(stm_threshold_10_90, "stm1090")
subset_weights_stm1585 = align_gene_weights(stm_threshold_15_85, "stm1585")
subset_weights_stm2575 = align_gene_weights(stm_threshold_25_75, "stm2575")

subset_weights_mixedstm1090 = align_gene_weights(mixed_stm_threshold_10_90, "mixed_stm1090")
subset_weights_mixedstm1585 = align_gene_weights(mixed_stm_threshold_15_85, "mixed_stm1585")
subset_weights_mixedstm2575 = align_gene_weights(mixed_stm_threshold_25_75, "mixed_stm2575")


# QC: Confirm all aligned gene-weight vectors have the same length.
subset_gene_weights = [
    subset_weights_stm1090,
    subset_weights_stm1585,
    subset_weights_stm2575,
    subset_weights_mixedstm1090,
    subset_weights_mixedstm1585,
    subset_weights_mixedstm2575,
]

gene_counts = [len(weights) for weights in subset_gene_weights]
if len(set(gene_counts)) == 1:
    print(f"All aligned gene-weight vectors have {gene_counts[0]} genes")
else:
    print("Gene counts differ between aligned gene-weight vectors")
    print(gene_counts)


# Convert gene weights to reaction weights with iMAT GPR parsing.
rxn_weights_stm1090 = gene_to_rxn_weights(control_model, subset_weights_stm1090)
rxn_weights_stm1585 = gene_to_rxn_weights(control_model, subset_weights_stm1585)
rxn_weights_stm2575 = gene_to_rxn_weights(control_model, subset_weights_stm2575)

rxn_weights_mixedstm1090 = gene_to_rxn_weights(control_model, subset_weights_mixedstm1090)
rxn_weights_mixedstm1585 = gene_to_rxn_weights(control_model, subset_weights_mixedstm1585)
rxn_weights_mixedstm2575 = gene_to_rxn_weights(control_model, subset_weights_mixedstm2575)


# QC: Confirm all reaction-weight vectors have the same length.
reaction_weights = [
    rxn_weights_stm1090,
    rxn_weights_stm1585,
    rxn_weights_stm2575,
    rxn_weights_mixedstm1090,
    rxn_weights_mixedstm1585,
    rxn_weights_mixedstm2575,
]

reaction_counts = [len(weights) for weights in reaction_weights]
if len(set(reaction_counts)) == 1:
    print(f"All reaction-weight vectors have {reaction_counts[0]} reactions")
else:
    print("Reaction counts differ between reaction-weight vectors")
    print(reaction_counts)


# Generate iMAT-constrained models from the control model and reaction weights.
stm1090_model = generate_model(
    model=control_model,
    rxn_weights=rxn_weights_stm1090,
    method="imat_restrictions",
)
stm1585_model = generate_model(
    model=control_model,
    rxn_weights=rxn_weights_stm1585,
    method="imat_restrictions",
)
stm2575_model = generate_model(
    model=control_model,
    rxn_weights=rxn_weights_stm2575,
    method="imat_restrictions",
)

mixed_stm1090_model = generate_model(
    model=control_model,
    rxn_weights=rxn_weights_mixedstm1090,
    method="imat_restrictions",
)
mixed_stm1585_model = generate_model(
    model=control_model,
    rxn_weights=rxn_weights_mixedstm1585,
    method="imat_restrictions",
)
mixed_stm2575_model = generate_model(
    model=control_model,
    rxn_weights=rxn_weights_mixedstm2575,
    method="imat_restrictions",
)


# QC: Print default medium exchange reactions, metabolites, and uptake values.
print("Default control-model medium:")
for rxn_id, uptake_value in control_model.medium.items():
    rxn = control_model.reactions.get_by_id(rxn_id)
    met = list(rxn.metabolites)[0]
    print(rxn_id, "|", met.id, "|", met.name, "|", uptake_value)


# -----------------------------------------------------------------------------
# Media definitions
# -----------------------------------------------------------------------------

# Base medium from the STM_V1_0 default medium.
default_medium = {
    "EX_ca2_e": 1000.0,
    "EX_cbl1_e": 0.01,
    "EX_cl_e": 1000.0,
    "EX_co2_e": 1000.0,
    "EX_cobalt2_e": 1000.0,
    "EX_cu2_e": 1000.0,
    "EX_mg2_e": 1000.0,
    "EX_fe2_e": 1000.0,
    "EX_fe3_e": 1000.0,
    "EX_mn2_e": 1000.0,
    "EX_mobd_e": 1000.0,
    "EX_na1_e": 1000.0,
    "EX_nh4_e": 1000.0,
    "EX_o2_e": 18.5,
    "EX_pi_e": 1000.0,
    "EX_zn2_e": 1000.0,
    "EX_so4_e": 1000.0,
    "EX_glc__D_e": 5.0,
    "EX_tungs_e": 1000.0,
    "EX_h2o_e": 100.0,
    "EX_h_e": 100.0,
    "EX_k_e": 1000.0,
}

# Infection medium: default medium plus infection-associated metabolites.
infected_medium = dict(default_medium)
infected_medium.update(
    {
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
        "EX_glc__D_e": 0.0,  # Remove glucose from the default medium.
    }
)

# Uninfected medium: lower tetrathionate and nitrate than infection medium.
uninfected_medium = dict(default_medium)
uninfected_medium.update(
    {
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
        "EX_glc__D_e": 0.0,  # Remove glucose from the default medium.
    }
)

base_media = {
    "default": default_medium,
    "infected": infected_medium,
    "uninfected": uninfected_medium,
}

original_models = {
    "control": control_model,
    "stm1090": stm1090_model,
    "stm1585": stm1585_model,
    "stm2575": stm2575_model,
    "mixed_stm1090": mixed_stm1090_model,
    "mixed_stm1585": mixed_stm1585_model,
    "mixed_stm2575": mixed_stm2575_model,
}

# Copy models before FBA so the original generated models are not altered.
fba_models = {model_name: model.copy() for model_name, model in original_models.items()}


# -----------------------------------------------------------------------------
# Oxygen sweep FBA
# -----------------------------------------------------------------------------

def run_oxygen_sweep(
    model_name: str,
    medium_condition: str,
    step_level: float = 0.5,
    o2_min: float = 0.0,
    o2_max: float = 18.5,
    o2_rxn: str = "EX_o2_e",
) -> pd.DataFrame:
    """Run FBA across an oxygen gradient for one model and one medium."""

    # Check that the requested medium and model exist.
    if medium_condition not in base_media:
        raise ValueError(
            f"{medium_condition} is not in base_media. "
            f"Available media: {list(base_media.keys())}"
        )

    if model_name not in fba_models:
        raise ValueError(
            f"{model_name} is not in fba_models. "
            f"Available models: {list(fba_models.keys())}"
        )

    if step_level <= 0:
        raise ValueError("step_level must be greater than 0.")

    # Build oxygen levels and include o2_max even if the step does not land on it.
    oxygen_levels = np.arange(o2_min, o2_max + step_level, step_level)
    oxygen_levels = oxygen_levels[oxygen_levels <= o2_max]
    oxygen_levels = np.round(oxygen_levels, 6)

    if oxygen_levels[-1] != o2_max:
        oxygen_levels = np.append(oxygen_levels, o2_max)

    base_medium = base_media[medium_condition]
    model = fba_models[model_name].copy()
    results = []

    # Run FBA once for each oxygen uptake value.
    for o2_value in oxygen_levels:
        medium = base_medium.copy()
        medium[o2_rxn] = float(o2_value)

        try:
            with model:
                model.medium = medium
                solution = model.optimize()

                objective_value = (
                    solution.objective_value if solution.status == "optimal" else np.nan
                )

                results.append(
                    {
                        "model": model_name,
                        "medium_condition": medium_condition,
                        "oxygen": float(o2_value),
                        "status": solution.status,
                        "objective_value": objective_value,
                    }
                )

        except Exception as error:
            results.append(
                {
                    "model": model_name,
                    "medium_condition": medium_condition,
                    "oxygen": float(o2_value),
                    "status": "failed",
                    "objective_value": np.nan,
                    "error": str(error),
                }
            )

    return pd.DataFrame(results).sort_values("oxygen").reset_index(drop=True)


def run_model_group(model_names: list[str], sample_type: str) -> pd.DataFrame:
    """Run oxygen sweeps for one sample group across all media."""

    group_results = []

    for model_name in model_names:
        for medium_condition in base_media:
            print(f"Running {model_name} in {medium_condition} medium")

            try:
                results_df = run_oxygen_sweep(
                    model_name=model_name,
                    medium_condition=medium_condition,
                    step_level=0.5,
                )
                results_df["sample_type"] = sample_type
                group_results.append(results_df)

            except Exception as error:
                print(f"FAILED: {model_name} in {medium_condition}")
                print(error)

                group_results.append(
                    pd.DataFrame(
                        [
                            {
                                "sample_type": sample_type,
                                "model": model_name,
                                "medium_condition": medium_condition,
                                "oxygen": np.nan,
                                "status": "failed_chunk",
                                "objective_value": np.nan,
                                "error": str(error),
                            }
                        ]
                    )
                )

    return pd.concat(group_results, ignore_index=True)


# Run the full oxygen sweep experiment.
control_oxygen_sweep_df = run_model_group(["control"], "control")
stm_oxygen_sweep_df = run_model_group(["stm1090", "stm1585", "stm2575"], "stm")
mixed_stm_oxygen_sweep_df = run_model_group(
    ["mixed_stm1090", "mixed_stm1585", "mixed_stm2575"],
    "mixed_stm",
)

all_oxygen_sweep_df = pd.concat(
    [control_oxygen_sweep_df, stm_oxygen_sweep_df, mixed_stm_oxygen_sweep_df],
    ignore_index=True,
)


# Save all oxygen sweep result tables to the shared output directory.
control_oxygen_sweep_df.to_csv(output_dir / "control_oxygen_sweep_df.csv", index=False)
stm_oxygen_sweep_df.to_csv(output_dir / "stm_oxygen_sweep_df.csv", index=False)
mixed_stm_oxygen_sweep_df.to_csv(
    output_dir / "mixed_stm_oxygen_sweep_df.csv",
    index=False,
)
all_oxygen_sweep_df.to_csv(output_dir / "all_oxygen_sweep_df.csv", index=False)

print(f"Saved CSV results to: {output_dir}")


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def save_oxygen_sweep_plot(
    df: pd.DataFrame,
    filter_column: str,
    filter_value: str,
    group_column: str,
    title: str,
    filename: str,
    xlim: tuple[float, float] = (0, 18.5),
    ylim: tuple[float, float] | None = (0, 1.8),
) -> None:
    """Plot and save oxygen sweep results for selected rows."""

    plot_df = df[df[filter_column] == filter_value].copy()

    # Remove full failed_chunk rows with no oxygen value.
    plot_df = plot_df.dropna(subset=["oxygen"])

    # Treat failed or missing objective values as no growth for plotting.
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
            label=group_name,
        )

        # Mark failed, missing, or zero-growth points in red.
        no_growth_df = group_df[
            (group_df["objective_value"].isna())
            | (np.isclose(group_df["plot_objective_value"], 0.0))
            | (group_df["status"] != "optimal")
        ]

        plt.scatter(
            no_growth_df["oxygen"],
            no_growth_df["plot_objective_value"],
            color="red",
            s=40,
            zorder=3,
        )

    plt.xlim(xlim)

    if ylim is not None:
        plt.ylim(ylim)

    plt.xlabel("Oxygen uptake value for EX_o2_e")
    plt.ylabel("Objective value")
    plt.title(title)
    plt.legend()
    plt.grid(True)

    save_path = output_dir / filename
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved plot: {save_path}")


# Plot the control model across all three media.
save_oxygen_sweep_plot(
    df=all_oxygen_sweep_df,
    filter_column="model",
    filter_value="control",
    group_column="medium_condition",
    title="Control model growth across oxygen gradient",
    filename="plot_01_control_three_media.png",
)

# Plot STM 10/90 across all three media.
save_oxygen_sweep_plot(
    df=all_oxygen_sweep_df,
    filter_column="model",
    filter_value="stm1090",
    group_column="medium_condition",
    title="STM 10/90 model growth across oxygen gradient",
    filename="plot_02_stm1090_three_media.png",
)

# Plot mixed STM 10/90 across all three media.
save_oxygen_sweep_plot(
    df=all_oxygen_sweep_df,
    filter_column="model",
    filter_value="mixed_stm1090",
    group_column="medium_condition",
    title="Mixed STM 10/90 model growth across oxygen gradient",
    filename="plot_03_mixed_stm1090_three_media.png",
)

# Compare control, STM 10/90, and mixed STM 10/90 in each medium.
selected_models = ["control", "stm1090", "mixed_stm1090"]
selected_model_df = all_oxygen_sweep_df[
    all_oxygen_sweep_df["model"].isin(selected_models)
].copy()

save_oxygen_sweep_plot(
    df=selected_model_df,
    filter_column="medium_condition",
    filter_value="default",
    group_column="model",
    title="Control, STM 10/90, and mixed STM 10/90 growth in default medium",
    filename="plot_04_three_models_default_medium.png",
)

save_oxygen_sweep_plot(
    df=selected_model_df,
    filter_column="medium_condition",
    filter_value="infected",
    group_column="model",
    title="Control, STM 10/90, and mixed STM 10/90 growth in infected medium",
    filename="plot_05_three_models_infected_medium.png",
)

save_oxygen_sweep_plot(
    df=selected_model_df,
    filter_column="medium_condition",
    filter_value="uninfected",
    group_column="model",
    title="Control, STM 10/90, and mixed STM 10/90 growth in uninfected medium",
    filename="plot_06_three_models_uninfected_medium.png",
)

print("Oxygen sweep script completed.")
