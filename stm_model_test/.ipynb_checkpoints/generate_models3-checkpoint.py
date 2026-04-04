import sys
import argparse
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
    p.add_argument("method")
    p.add_argument("--output", required = True)
    args = p.parse_args()

    # set solver 
    Configuration().solver = "gurobi"
    Configuration().solver_configuration = {"Threads": args.threads,}

    # load the test S.Tm LT2 GEM in cobrapy
    base_model = load_model("STM_v1_0")
    base_model.solver = "gurobi"

    # read in the csv files of the threshold rna-seq data
    threshold_csv = pd.read_csv(args.threshold_csv)

    # python set
    threshold_genes = set(threshold_csv["genes"].astype(str).str.strip())

    # make a list of genes and save it as a set
    model_genes = set([g.id for g in base_model.genes])

    # find genes that exist in the model and RNA-seq
    overlap = threshold_genes.intersection(model_genes)

    # find genes that appear in the RNA-Seq data but are not in the model
    only_stm = threshold_genes.difference(model_genes)

    # find genes that exist in the model but do not appear in RNA-Seq data
    only_model = model_genes.difference(threshold_genes)

    # subset the rna-seq expression data to only have the genes present in the model
    threshold_model = threshold_csv[threshold_csv["genes"].isin(model_genes)]

    # to this subset - add the 20 genes present in the model but missing in the rna-seq data and assign a value of 0
    # convert from dataframe to panda series, make sure all the genes from the model exist, and fill missing values with 0
    expression_subset = threshold_model.set_index("genes")["threshold_category"].reindex(model_genes, fill_value=0)
    
    # convert gene weights to reaction weights
    rxn_weights = gene_to_rxn_weights(base_model, expression_subset)

    try:
        print(f"{args.method} generating model for {args.threshold_csv}")
        model = generate_model(model = base_model, rxn_weights = rxn_weights, method = args.method)
        save_json_model(
            model,
            args.output
            )
    except Exception as e:
        print(f"FAILED SUBSET for {args.threshold_csv}: {e}")

if __name__ == '__main__':
    sys.exit(main())