# core
import numpy as np
import pandas as pd

# COBRApy
from cobra.core.configuration import Configuration
from cobra.io import load_model
from cobra.io import save_json_model

#iMATpy
from imatpy.imat import imat
from imatpy.model_creation import generate_model
from imatpy.parse_gpr import gene_to_rxn_weights

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-t', '--threads', default=8, type = int,
                   help='number of threads to use')
    p.add_argument("threshold_csv")
    args = p.parse_args()


    # set solver 
    Configuration().solver = "gurobi"
    Configuration().solver_configuration = {"Threads": args.threads,}

    # load the test S.Tm LT2 GEM in cobrapy
    base_model = load_model("STM_v1_0")
    base_model.solver = "gurobi"

    # read in the csv files of the threshold rna-seq data
    stm_threshold_10_90 = pd.read_csv(args.threshold_csv)

    # python set
    threshold_genes = set(stm_threshold_10_90["genes"].astype(str).str.strip())

    # make a list of genes and save it as a set
    model_genes = set([g.id for g in base_model.genes])

    # find genes that exist in the model and RNA-seq
    overlap = threshold_genes.intersection(model_genes)

    # find genes that appear in the RNA-Seq data but are not in the model
    only_stm = threshold_genes.difference(model_genes)

    # find genes that exist in the model but do not appear in RNA-Seq data
    only_model = model_genes.difference(threshold_genes)

    # subset the rna-seq expression data to only have the genes present in the model
    stm_threshold_10_90_model = stm_threshold_10_90[stm_threshold_10_90["genes"].isin(model_genes)]
    stm_threshold_15_85_model = stm_threshold_15_85[stm_threshold_15_85["genes"].isin(model_genes)]
    stm_threshold_25_75_model = stm_threshold_25_75[stm_threshold_25_75["genes"].isin(model_genes)]

    mixed_stm_threshold_10_90_model = mixed_stm_threshold_10_90[mixed_stm_threshold_10_90["genes"].isin(model_genes)]
    mixed_stm_threshold_15_85_model = mixed_stm_threshold_15_85[mixed_stm_threshold_15_85["genes"].isin(model_genes)]
    mixed_stm_threshold_25_75_model = mixed_stm_threshold_25_75[mixed_stm_threshold_25_75["genes"].isin(model_genes)]

    # to these subsets add the 20 genes present in the model but missing in the rna-seq data and assign a value of 0
    # convert from dataframe to panda series, make sure all the genes from the model exist, and fill missing values with 0

    stm_10_90_expr = stm_threshold_10_90_model.set_index("genes")["threshold_category"].reindex(model_genes, fill_value=0)
    stm_15_85_expr = stm_threshold_15_85_model.set_index("genes")["threshold_category"].reindex(model_genes, fill_value=0)
    stm_25_75_expr = stm_threshold_25_75_model.set_index("genes")["threshold_category"].reindex(model_genes, fill_value=0)

    mixed_10_90_expr = mixed_stm_threshold_10_90_model.set_index("genes")["threshold_category"].reindex(model_genes, fill_value=0)
    mixed_15_85_expr = mixed_stm_threshold_15_85_model.set_index("genes")["threshold_category"].reindex(model_genes, fill_value=0)
    mixed_25_75_expr = mixed_stm_threshold_25_75_model.set_index("genes")["threshold_category"].reindex(model_genes, fill_value=0)

    # convert gene weights to reaction weights
    stm_10_90_rxn = gene_to_rxn_weights(base_model, stm_10_90_expr)
    stm_15_85_rxn = gene_to_rxn_weights(base_model, stm_15_85_expr)
    stm_25_75_rxn = gene_to_rxn_weights(base_model, stm_25_75_expr)

    mixed_10_90_rxn = gene_to_rxn_weights(base_model, mixed_10_90_expr)
    mixed_15_85_rxn = gene_to_rxn_weights(base_model, mixed_15_85_expr)
    mixed_25_75_rxn = gene_to_rxn_weights(base_model, mixed_25_75_expr)

    # make a dictionary of reaction weights
    reaction_weights = {
        "stm_10_90": stm_10_90_rxn,
        "stm_15_85": stm_15_85_rxn,
        "stm_25_75": stm_25_75_rxn,
        "mixed_10_90": mixed_10_90_rxn,
        "mixed_15_85": mixed_15_85_rxn,
        "mixed_25_75": mixed_25_75_rxn
    }

    # Run iMAT 
    imat_results = {}

    for name, rxn_weights in reaction_weights.items():
        print(f"Running iMAT for {name}")
        imat_result = imat(model = base_model, rxn_weights = rxn_weights)
        imat_results[name] = imat_result
        
    # Save the flux distribution 
    for name, result in imat_results.items():
        flux_df = pd.DataFrame({
            "reaction": list(result.fluxes.index),
            "flux": list(result.fluxes.values)
        })
        flux_df.to_csv(
            f"/home/gmvaz/2026_GEMs/stm_model_test/run_031026/fluxes/{name}_imat_fluxes.csv",
            index=False
        )

    # Generate context models 
    # With imat_restrictions
    context_model_imatrest = {}

    for name, rxn_weights in reaction_weights.items():
        try:
            print(f"[imat_restrictions] generating model for {name}")
            model = generate_model(model = base_model, rxn_weights = rxn_weights, method = 'imat_restrictions')
            context_model_imatrest[name] = model
        except Exception as e:
            print(f"FAILED IMAT RESTRICTIONS for {name}: {e}")


    # With simple_bounds
    context_model_simpbounds = {}

    for name, rxn_weights in reaction_weights.items():
        try:
            print(f"[simple_bounds] generating model for {name}")
            model = generate_model(model = base_model, rxn_weights = rxn_weights, method = 'simple_bounds')
            context_model_simpbounds[name] = model
        except Exception as e:
            print(f"FAILED SIMPLE BOUNDS for {name}: {e}")

    # With eliminiate_below_threshold
    context_model_subset = {}

    for name, rxn_weights in reaction_weights.items():
        try:
            print(f"[eliminate_below_threshold] generating model for {name}")
            model = generate_model(model = base_model, rxn_weights = rxn_weights, method = 'eliminate_below_threshold')
            context_model_subset[name] = model
        except Exception as e:
            print(f"FAILED SUBSET for {name}: {e}")

    # With FVA
    context_model_FVA = {}

    for name, rxn_weights in reaction_weights.items():
        try:
            print(f"[FVA] generating model for {name}")
            model = generate_model(model = base_model, rxn_weights = rxn_weights, method = 'FVA')
            context_model_FVA[name] = model
        except Exception as e:
            print(f"FAILED FVA for {name}: {e}")

    # With MILP
    context_model_MILP = {}

    for name, rxn_weights in reaction_weights.items():
        try:
            print(f"[MILP] generating model for {name}")
            model = generate_model(base_model, rxn_weights, method='MILP')
            context_model_MILP[name] = model
        except Exception as e:
            print(f"FAILED MILP for {name}: {e}")


    # Make a dictionary of the model dictionaries
    model_sets = {
        "imatrest": context_model_imatrest,
        "simple": context_model_simpbounds,
        "subset": context_model_subset,
        "FVA": context_model_FVA,
        "MILP": context_model_MILP
    }

    # Save the models 
    for method, models in model_sets.items():
        for name, model in models.items():
            save_json_model(
                model,
                f"/home/gmvaz/2026_GEMs/stm_model_test/run_031026/context_models/{name}_{method}.json"
            )

if __name__ == '__main__':
    sys.exit(main())