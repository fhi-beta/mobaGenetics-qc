## GWAS test pipeline

This pipeline conducts a test GWAS on the MoBa genotypes.

## Content

- [docs](docs): Documentation produced by the pipeline on the different releases
- [envs](envs): Files used to manage conda environments
- [utils](utils): Utilities scripts
- [analysis.yaml](analysis.yaml): Parameters for the analysis
- [gwas.snake](gwas.snake): Main pipeline file
- [prepare_phenotypes_file.R](prepare_phenotypes_file.R): Script used to prepare the phenotypes for regenie
- [association_qc_no_cojo.R](association_qc_no_cojo.R): Script used to write the documentation on the association results
