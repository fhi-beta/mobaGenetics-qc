# Fam file reconstruction in snp017d
- Number of samples in the genotyping data: 5573.
## Samples not in Medical Birth Regsitry
11 samples with missing birth year, assumed to be parent in the following.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 9 |
| Parent-offspring| 460 |
| Full siblings| 26 |
| 2nd degree| 0 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](fam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 24 |
| Male | 1 |
| Female | 2144 |

![](fam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 1110 |
| Female | 2 |

![](fam_reconstruction/father_sex_plot.png)
## Children sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 6 |
| Male | 1103 |
| Female | 1183 |

![](fam_reconstruction/children_sex_plot.png)
## Parental relationships
11 sentrix IDs missing from ID file. These are not counted as individuals.
###  Individuals
5562 individuals in total. Breakdown excluding multiple same-sex parents:
 -  395 children
 -  258 mothers
 -  199 fathers
 -  259 mother-child pairs
 -  200 father-child pairs
 -  64 mother-father-child trios
 -  4710 unrelated

Multiple same-sex parents (at the individual level):
 -  0 children with more than one mother detected
 -  0 children with more than one father detected
 -  0 children with more than one mother in registry
 -  0 children with more than one father in registry

261 mother-child relationships expected.
- 259 (99.23%) recovered by genetic relationships.
- 2 (0.77%) not recovered by genetic relationships.


199 father-child relationships expected.
- 199 (100%) recovered by genetic relationships.
- 0 (0%) not recovered by genetic relationships.


259 mother-child relationships detected.
- 259 (100%) matched to registry.
- 0 (0%) not matched to registry.


200 father-child relationships detected.
- 199 (99.5%) matched to registry.
- 1 (0.5%) not matched to registry.


###  Samples
5573 samples in total. Breakdown excluding multiple same-sex parents:
 -  395 children
 -  258 mothers
 -  199 fathers
 -  259 mother-child pairs
 -  200 father-child pairs
 -  64 mother-father-child trios
 -  4721 unrelated

Multiple same-sex parents (at the sample level):
 -  0 children with more than one mother detected
 -  0 children with more than one father detected
 -  3883 children with more than one mother in registry
 -  2490 children with more than one father in registry

261 mother-child relationships expected.
- 259 (99.23%) recovered by genetic relationships.
- 2 (0.77%) not recovered by genetic relationships.


199 father-child relationships expected.
- 199 (100%) recovered by genetic relationships.
- 0 (0%) not recovered by genetic relationships.


259 mother-child relationships detected.
- 259 (100%) matched to registry.
- 0 (0%) not matched to registry.


200 father-child relationships detected.
- 199 (99.5%) matched to registry.
- 1 (0.5%) not matched to registry.


## Exclusion
- Number of samples excluded: 7
