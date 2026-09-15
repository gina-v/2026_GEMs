MODEL = ["stm_1090", "stm_1585", "stm_2575", "mixed_stm1090", "mixed_stm1585", "mixed_stm2575"]

MEDIUM = ["A_infection_medium", "B_infection_medium", "C_infection_medium", "D_infection_medium", "E_infection_medium", "F_infection_medium", "G_infection_medium"]

RXN = ["EX_tet_e", "EX_no3_e", "EX_etha_e", "EX_12ppd__R_e", "EX_12ppd__S_e", "EX_lac__L_e", "EX_for_e", "EX_succ_e", "EX_asp__L_e", "EX_mal__L_e", "EX_galct__D_e", "EX_glcr_e", "EX_pyr_e", "EX_arab__L_e", "EX_but_e", "EX_ac_e", "EX_ppa_e"]

# Global variables for metabolite_sweep.py
MODEL2 = ["stm_1090", "stm_1585", "stm_2575", "mixed_stm1090", "mixed_stm1585", "mixed_stm2575"]

MEDIUM2 = ["A_infection_medium", "B_infection_medium", "D_infection_medium", "F_infection_medium"]

RXN2 = ['EX_alaala_e', 'EX_all__D_e', 'EX_alltn_e', 'EX_amp_e', 'EX_anhgm_e']

rule all:
    input:
        expand("nutrient_sweep/{model}.{medium}.{rxn}.csv", model = MODEL, medium = MEDIUM, rxn = RXN),
        expand("metabolite_sweep/{model}.{medium}.{rxn}.csv", model = MODEL2, medium = MEDIUM2, rxn = RXN2)

rule make_csv:
    output:
        "nutrient_sweep/{model}.{medium}.{rxn}.csv"
    shell:
        "python nutrient_sweep.py {wildcards.model} {wildcards.medium} {wildcards.rxn} {output}"

rule make_csv2:
    output:
        "metabolite_sweep/{model}.{medium}.{rxn}.csv"
    shell:
        "python metabolite_sweep.py {wildcards.model} {wildcards.medium} {wildcards.rxn} {output}"

# Notes from 09/01/26 meeting with TB
# Copy nutrient_sweep.py to metabolite_sweep.py, edit the sweep range, 
# RXN will only have reactions that are in the infection medium 
# RNX2 will have all other metabolites from STM_V1_0
# 

# snakemake -s Snakefile.py -j Cores_per_job -k
# snakemake -k will keep running even if there's an error

# MODEL LIST: ["stm_1090", "stm_1585", "stm_2575", "mixed_stm1090", "mixed_stm1090", "mixed_stm1585", "mixed_stm2575"]

# MEDIUM LIST: MEDIUM = ["A_infection_medium", "B_infection_medium", "C_infection_medium", "D_infection_medium", "E_infection_medium", "F_infection_medium", "G_infection_medium"]

# RXN LIST: ["EX_tet_e", "EX_no3_e", "EX_etha_e", "EX_12ppd__R_e", "EX_12ppd__S_e", "EX_lac__L_e", "EX_for_e", "EX_succ_e", "EX_asp__L_e", "EX_mal__L_e", "EX_galct__D_e", "EX_glcr_e", "EX_pyr_e", "EX_arab__L_e", "EX_but_e", "EX_ac_e", "EX_ppa_e"]