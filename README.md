# Strong isolation by distance and evidence of population microstructure reflect ongoing Plasmodium falciparum transmission in Zanzibar

This repository is for the manuscript entitled ["Strong isolation by distance and evidence of population microstructure reflect ongoing Plasmodium falciparum transmission in Zanzibar"](https://elifesciences.org/reviewed-preprints/90173).

The MIPAnalyzer objects and vcfR objects are in the [data/data_objects](https://github.com/sconnelly007/TAN_MIP/tree/main/data/data_objects) directory. These objects are organized as follows:
- `biallelic_distances.rds`: the MIPAnalyzer object from Verity et al. 2020
- `biallelic_processed_Kenya.rds`: the MIPAnalyzer object from our dataset from Kenya generated in this study.
- `biallelic_processed_TAN.rds `: the MIPAnalyzer object from our dataset from Zanzibar generated in this study. Metadata is in the [metadata](https://github.com/sconnelly007/TAN_MIP/tree/main/data/metadata) directory.
- `vcf_nomiss.vcf.gz`: the vcfR object from our dataset from Zanzibar generated in this study. This vcf was filtered to targeted sites, by sample and loci missingness, and on minor and major allele frequency. Metadata is in the [metadata](https://github.com/sconnelly007/TAN_MIP/tree/main/data/metadata) directory.


The R scripts are in the [R](https://github.com/sconnelly007/TAN_MIP/tree/main/R) directory. These scripts are organized as follows:
- IBD
    - `00_IBD_rev.Rmd`: This script is used to calulate IBD and perform downsteam analysis, such as PCA, DPAC, network analysis and isolation by distance.
    - `01_PCA_Combined_rev.Rmd`: This script is used to perform PCA analysis on malaria genetic data from throughout Africa to put the samples in this study in context.
    - `02_COI_Fws.Rmd`: This script calculated COI and Fws from the study's data object, which is saved [here](https://github.com/sconnelly007/TAN_MIP/tree/main/data/data_objects/vcf_nomiss.vcf.gz). The output from the THEREALMcCOIL is saved in [data/COI_output](https://github.com/sconnelly007/TAN_MIP/tree/main/data/COI_output).
- DR2_Analysis
    - `00_prevalence_table.R`: Using the genotypes data table in the [tables](https://github.com/sconnelly007/TAN_MIP/tree/main/data) directory, this script processes and creates tables for the prevalence of different polymorphisms related to malarial drug resistance. 

