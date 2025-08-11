## Test GWAS

This release was tested by running a test GWAS of the children's birth weight against the genotypes of chidren, mothers, and fathers separately.

### Results

A GWAS of birth weight has been run for all children using a [test GWAS pipeline](../../readme.md). Briefly, extreme outliers were excluded, and a GWAS was conducted using [regenie](rgcgithub.github.io/regenie). Sex and genotyping batch, and ten PCs were used as default covariates.

- [Children](pc_covar/regenie_no_cojo/weight_birth/pop_children_pheno_birth_weight.md): GWAS of the children birth weight against the children's genome.
- [Mothers](pc_covar/regenie_no_cojo/weight_birth/pop_mothers_pheno_birth_weight.md): GWAS of the children birth weight against the mothers' genome.
- [Mothers](pc_covar/regenie_no_cojo/weight_birth/pop_fathers_pheno_birth_weight.md): GWAS of the children birth weight against the fathers' genome.

### Disclaimer

This GWAS is run for testing puroses only. The results should not be used for other applications. Note that phenotypes were not curated and analyses do not account for covariates like gestational duration, hence limiting statistical power.

