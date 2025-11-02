# Fam file reconstruction in snp011
- Number of samples in the genotyping data: 5182.
## Samples not in Medical Birth Regsitry
37 samples with missing birth year, assumed to be parent in the following.
## Relationship inference
| Relationship |   |
| ------------ | - |
| Duplicates or monozygotic twins| 5 |
| Parent-offspring| 2028 |
| Full siblings| 28 |
| 2nd degree| 0 |
| 3rd degree| 0 |
| 4th degree| 0 |
| Unrelated| 0 |

![](fam_reconstruction/ibd_plot.png)
## Mother sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 34 |
| Male | 1 |
| Female | 1430 |

![](fam_reconstruction/mother_sex_plot.png)
## Father sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 1 |
| Male | 772 |
| Female | 2 |

![](fam_reconstruction/father_sex_plot.png)
## Children sex check
| Inferred sex |   |
| ------------ | - |
| Unknown | 13 |
| Male | 2114 |
| Female | 815 |

![](fam_reconstruction/children_sex_plot.png)
## Parental relationships
37 sentrix IDs missing from ID file. These are not counted as individuals.
###  Individuals
5145 individuals in total. Breakdown excluding multiple same-sex parents:
 -  1378 children
 -  1312 mothers
 -  689 fathers
 -  1328 mother-child pairs
 -  698 father-child pairs
 -  648 mother-father-child trios
 -  1766 unrelated

Multiple same-sex parents (at the individual level):
 -  0 children with more than one mother detected
 -  0 children with more than one father detected
 -  0 children with more than one mother in registry
 -  0 children with more than one father in registry

1335 mother-child relationships expected.
- 1326 (99.33%) recovered by genetic relationships.
- 9 (0.67%) not recovered by genetic relationships.


703 father-child relationships expected.
- 691 (98.29%) recovered by genetic relationships.
- 12 (1.71%) not recovered by genetic relationships.


1328 mother-child relationships detected.
- 1326 (99.85%) matched to registry.
- 2 (0.15%) not matched to registry.


698 father-child relationships detected.
- 691 (99%) matched to registry.
- 7 (1%) not matched to registry.


###  Samples
5182 samples in total. Breakdown excluding multiple same-sex parents:
 -  1378 children
 -  1312 mothers
 -  689 fathers
 -  1328 mother-child pairs
 -  698 father-child pairs
 -  648 mother-father-child trios
 -  1803 unrelated

Multiple same-sex parents (at the sample level):
 -  0 children with more than one mother detected
 -  0 children with more than one father detected
 -  3883 children with more than one mother in registry
 -  2490 children with more than one father in registry

1335 mother-child relationships expected.
- 1326 (99.33%) recovered by genetic relationships.
- 9 (0.67%) not recovered by genetic relationships.


703 father-child relationships expected.
- 691 (98.29%) recovered by genetic relationships.
- 12 (1.71%) not recovered by genetic relationships.


1328 mother-child relationships detected.
- 1326 (99.85%) matched to registry.
- 2 (0.15%) not matched to registry.


698 father-child relationships detected.
- 691 (99%) matched to registry.
- 7 (1%) not matched to registry.


## Exclusion
- Number of samples excluded: 27
