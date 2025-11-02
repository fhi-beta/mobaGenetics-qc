# Fam file reconstruction in snp001
- Number of samples in the genotyping data: 18492.
## Samples not in Medical Birth Regsitry
140 samples with missing birth year, assumed to be parent in the following.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 142 |
| Parent-offspring| 11223 |
| Full siblings| 196 |
| 2nd degree| 0 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](fam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 20 |
| Male | 2 |
| Female | 6193 |

![](fam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 0 |
| Male | 6150 |
| Female | 6 |

![](fam_reconstruction/father_sex_plot.png)
## Children sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 3 |
| Male | 3143 |
| Female | 2975 |

![](fam_reconstruction/children_sex_plot.png)
## Parental relationships
140 sentrix IDs missing from ID file. These are not counted as individuals.
###  Individuals
18220 individuals in total. Breakdown excluding multiple same-sex parents:
 -  5704 children
 -  5481 mothers
 -  5375 fathers
 -  5530 mother-child pairs
 -  5425 father-child pairs
 -  5251 mother-father-child trios
 -  1660 unrelated

Multiple same-sex parents (at the individual level):
 -  5 children with more than one mother detected
 -  5 children with more than one father detected
 -  0 children with more than one mother in registry
 -  0 children with more than one father in registry

5532 mother-child relationships expected.
- 5530 (99.96%) recovered by genetic relationships.
- 2 (0.04%) not recovered by genetic relationships.


5431 father-child relationships expected.
- 5425 (99.89%) recovered by genetic relationships.
- 6 (0.11%) not recovered by genetic relationships.


5535 mother-child relationships detected.
- 5530 (99.91%) matched to registry.
- 5 (0.09%) not matched to registry.


5430 father-child relationships detected.
- 5425 (99.91%) matched to registry.
- 5 (0.09%) not matched to registry.


###  Samples
18492 samples in total. Breakdown excluding multiple same-sex parents:
 -  5774 children
 -  5513 mothers
 -  5400 fathers
 -  5603 mother-child pairs
 -  5494 father-child pairs
 -  5323 mother-father-child trios
 -  1807 unrelated

Multiple same-sex parents (at the sample level):
 -  61 children with more than one mother detected
 -  56 children with more than one father detected
 -  3883 children with more than one mother in registry
 -  2490 children with more than one father in registry

5623 mother-child relationships expected.
- 5621 (99.96%) recovered by genetic relationships.
- 2 (0.04%) not recovered by genetic relationships.


5518 father-child relationships expected.
- 5512 (99.89%) recovered by genetic relationships.
- 6 (0.11%) not recovered by genetic relationships.


5664 mother-child relationships detected.
- 5621 (99.24%) matched to registry.
- 43 (0.76%) not matched to registry.


5550 father-child relationships detected.
- 5512 (99.32%) matched to registry.
- 38 (0.68%) not matched to registry.


## Exclusion
- Number of samples excluded: 19
