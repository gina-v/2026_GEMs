MODEL = ["stm_1090", "stm_1585", "mixed_stm1090", "mixed_stm1585"]

MEDIUM = ["A_infection_medium", "D_infection_medium", "F_infection_medium"]

RXN = ["EX_tet_e", "EX_no3_e", "EX_etha_e"]

rule all:
    input:
        expand("nutrient_sweep/{model}.{medium}.{rxn}.csv", model = MODEL, medium = MEDIUM, rxn = RXN)

rule make_csv:
    output:
        "nutrient_sweep/{model}.{medium}.{rxn}.csv"
    shell:
        "python nutrient_sweep.py {wildcards.model} {wildcards.medium} {wildcards.rxn} {output}"




# snakemake -s Snakefile.py -j Cores_per_job -k
# -k will keep running even if there's an error

# MODEL LIST: ["stm_1090", "stm_1585", "stm_2575", "mixed_stm1090", "mixed_stm1090", "mixed_stm1585", "mixed_stm2575"]

# MEDIUM LIST: MEDIUM = ["A_infection_medium", "B_infection_medium", "C_infection_medium", "D_infection_medium", "E_infection_medium", "F_infection_medium", "G_infection_medium"]

# RXN LIST: ["EX_tet_e", "EX_no3_e", "EX_etha_e", "EX_12ppd__R_e", "EX_12ppd__S_e", "EX_lac__L_e", "EX_for_e", "EX_succ_e", "EX_asp__L_e", "EX_mal__L_e", "EX_galct__D_e", "EX_glcr_e", "EX_pyr_e", "EX_arab__L_e", "EX_but_e", "EX_ac_e", "EX_ppa_e"]