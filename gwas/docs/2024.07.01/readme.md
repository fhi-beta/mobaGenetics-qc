## Test GWAS

The present release was tested by running a test GWAS against birth weight.

### Disclaimer

This GWAS is run for testing puroses only. The results should not be used for other purposes. Note that phenotypes were not curated and analyses do not account for covariates like gestational duration, hence limiting statistical power.

### Results

A GWAS of birth weight has been run for all children using a [test GWAS pipeline](../../gwas.snake). Briefly, extreme outliers were excluded, and a GWAS was conducted using [regenie](rgcgithub.github.io/regenie). Sex and genotyping batch were used as default covariates.

- [All samples](regenie_no_cojo/weight_birth/pop_children_pheno_birth_weight.md): GWAS of the birth weight for all children.
- [All samples with PCs](pc_covar/regenie_no_cojo/weight_birth/pop_children_pheno_birth_weight.md): GWAS of the birth weight for all children using 10 PCs as covariates.

