# Fam file reconstruction in snp002
- Number of samples in the genotyping data: 1665.
## Samples not in Medical Birth Regsitry
5 samples with missing birth year, assumed to be parent in the following.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 10 |
| Parent-offspring| 262 |
| Full siblings| 3 |
| 2nd degree| 0 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](fam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 2 |
| Female | 468 |

![](fam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 462 |
| Female | 0 |

![](fam_reconstruction/father_sex_plot.png)
## Children sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 400 |
| Female | 333 |

![](fam_reconstruction/children_sex_plot.png)
## Parental relationships
5 sentrix IDs missing from ID file. These are not counted as individuals.
###  Individuals
1650 individuals in total. Breakdown excluding multiple same-sex parents:
 -  179 children
 -  162 mothers
 -  96 fathers
 -  162 mother-child pairs
 -  96 father-child pairs
 -  79 mother-father-child trios
 -  1213 unrelated

Multiple same-sex parents (at the individual level):
 -  0 children with more than one mother detected
 -  0 children with more than one father detected
 -  0 children with more than one mother in registry
 -  0 children with more than one father in registry

163 mother-child relationships expected.
- 162 (99.39%) recovered by genetic relationships.
- 1 (0.61%) not recovered by genetic relationships.


96 father-child relationships expected.
- 96 (100%) recovered by genetic relationships.
- 0 (0%) not recovered by genetic relationships.


162 mother-child relationships detected.
- 162 (100%) matched to registry.
- 0 (0%) not matched to registry.


96 father-child relationships detected.
- 96 (100%) matched to registry.
- 0 (0%) not matched to registry.


###  Samples
1665 samples in total. Breakdown excluding multiple same-sex parents:
 -  181 children
 -  162 mothers
 -  96 fathers
 -  164 mother-child pairs
 -  97 father-child pairs
 -  80 mother-father-child trios
 -  1226 unrelated

Multiple same-sex parents (at the sample level):
 -  1 children with more than one mother detected
 -  0 children with more than one father detected
 -  3883 children with more than one mother in registry
 -  2490 children with more than one father in registry

166 mother-child relationships expected.
- 165 (99.4%) recovered by genetic relationships.
- 1 (0.6%) not recovered by genetic relationships.


97 father-child relationships expected.
- 97 (100%) recovered by genetic relationships.
- 0 (0%) not recovered by genetic relationships.


165 mother-child relationships detected.
- 165 (100%) matched to registry.
- 0 (0%) not matched to registry.


97 father-child relationships detected.
- 97 (100%) matched to registry.
- 0 (0%) not matched to registry.


## Exclusion
- Number of samples excluded: 3
