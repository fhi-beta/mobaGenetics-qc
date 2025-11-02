# Fam file reconstruction in snp017b
- Number of samples in the genotyping data: 4651.
## Samples not in Medical Birth Regsitry
21 samples with missing birth year, assumed to be parent in the following.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 3 |
| Parent-offspring| 219 |
| Full siblings| 22 |
| 2nd degree| 0 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](fam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 12 |
| Male | 2 |
| Female | 2177 |

![](fam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 1013 |
| Female | 1 |

![](fam_reconstruction/father_sex_plot.png)
## Children sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 1 |
| Male | 738 |
| Female | 707 |

![](fam_reconstruction/children_sex_plot.png)
## Parental relationships
21 sentrix IDs missing from ID file. These are not counted as individuals.
###  Individuals
4630 individuals in total. Breakdown excluding multiple same-sex parents:
 -  210 children
 -  193 mothers
 -  21 fathers
 -  198 mother-child pairs
 -  21 father-child pairs
 -  9 mother-father-child trios
 -  4206 unrelated

Multiple same-sex parents (at the individual level):
 -  0 children with more than one mother detected
 -  0 children with more than one father detected
 -  0 children with more than one mother in registry
 -  0 children with more than one father in registry

199 mother-child relationships expected.
- 197 (98.99%) recovered by genetic relationships.
- 2 (1.01%) not recovered by genetic relationships.


19 father-child relationships expected.
- 19 (100%) recovered by genetic relationships.
- 0 (0%) not recovered by genetic relationships.


198 mother-child relationships detected.
- 197 (99.49%) matched to registry.
- 1 (0.51%) not matched to registry.


21 father-child relationships detected.
- 19 (90.48%) matched to registry.
- 2 (9.52%) not matched to registry.


###  Samples
4651 samples in total. Breakdown excluding multiple same-sex parents:
 -  210 children
 -  193 mothers
 -  21 fathers
 -  198 mother-child pairs
 -  21 father-child pairs
 -  9 mother-father-child trios
 -  4227 unrelated

Multiple same-sex parents (at the sample level):
 -  0 children with more than one mother detected
 -  0 children with more than one father detected
 -  3883 children with more than one mother in registry
 -  2490 children with more than one father in registry

199 mother-child relationships expected.
- 197 (98.99%) recovered by genetic relationships.
- 2 (1.01%) not recovered by genetic relationships.


19 father-child relationships expected.
- 19 (100%) recovered by genetic relationships.
- 0 (0%) not recovered by genetic relationships.


198 mother-child relationships detected.
- 197 (99.49%) matched to registry.
- 1 (0.51%) not matched to registry.


21 father-child relationships detected.
- 19 (90.48%) matched to registry.
- 2 (9.52%) not matched to registry.


## Exclusion
- Number of samples excluded: 6
